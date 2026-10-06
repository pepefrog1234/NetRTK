/*------------------------------------------------------------------------------
 * gga.cpp : NMEA GGA parser and base64 helpers for the NTRIP front end
 *
 *   $xxGGA,hhmmss.ss,ddmm.mmmm,N,dddmm.mmmm,E,q,nn,hdop,alt,M,sep,M,age,refid*hh
 *   field   1         2        3 4          5 6 7  8    9   10 11 12 13  14
 *----------------------------------------------------------------------------*/
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include "nrtk/ntrip_front.hpp"

namespace nrtk {

namespace {

/* split at ',' keeping empty fields ------------------------------------------*/
std::vector<std::string> split_fields(const std::string &s)
{
    std::vector<std::string> v;
    for (size_t p=0;;) {
        size_t q=s.find(',',p);
        if (q==std::string::npos) {v.push_back(s.substr(p)); break;}
        v.push_back(s.substr(p,q-p));
        p=q+1;
    }
    return v;
}
int hexval(char c)
{
    if (c>='0'&&c<='9') return c-'0';
    if (c>='A'&&c<='F') return c-'A'+10;
    if (c>='a'&&c<='f') return c-'a'+10;
    return -1;
}
/* decimal number [+-]ddd[.ddd] (locale independent, no exponent) -------------*/
bool to_num(const std::string &s, double &v)
{
    const char *p=s.c_str();
    bool neg=false;
    double m=0.0;
    int nd=0,nf=0;

    if (*p=='+'||*p=='-') {neg=*p=='-'; p++;}
    for (;*p>='0'&&*p<='9';p++,nd++) m=m*10.0+(*p-'0');
    if (*p=='.') for (p++;*p>='0'&&*p<='9';p++,nf++) m=m*10.0+(*p-'0');
    if (*p||nd+nf==0||nd+nf>17) return false;
    v=(neg?-m:m)/std::pow(10.0,nf);
    return true;
}
bool to_int(const std::string &s, int &v)
{
    if (s.empty()||s.size()>6) return false;
    v=0;
    for (char c:s) {
        if (c<'0'||c>'9') return false;
        v=v*10+(c-'0');
    }
    return true;
}
/* ddmm.mmmm (lat) / dddmm.mmmm (lon) + hemisphere -> rad -----------------------*/
bool to_angle(const std::string &val, const std::string &hemi, double maxdeg,
              char pos, char neg, double &rad)
{
    double v,deg,min;

    if (val[0]<'0'||val[0]>'9'||!to_num(val,v)) return false;
    deg=std::floor(v/100.0);
    min=v-deg*100.0;
    if (min>=60.0||deg+min/60.0>maxdeg) return false;
    if (hemi.size()!=1||(hemi[0]!=pos&&hemi[0]!=neg)) return false;
    rad=(deg+min/60.0)*D2R*(hemi[0]==neg?-1.0:1.0);
    return true;
}
/* hhmmss[.sss] -> time of day (s) ---------------------------------------------*/
bool to_tod(const std::string &s, double &tod)
{
    int hh,mm;
    double ss;

    if (s.size()<6||s[4]<'0'||s[4]>'9'||s[5]<'0'||s[5]>'9') return false;
    if (!to_int(s.substr(0,2),hh)||!to_int(s.substr(2,2),mm)||!to_num(s.substr(4),ss)) return false;
    if (hh>23||mm>59||ss>=61.0) return false;
    tod=hh*3600.0+mm*60.0+ss;
    return true;
}
bool is_metre(const std::vector<std::string> &f, size_t i)
{
    return i>=f.size()||f[i].empty()||f[i]=="M"||f[i]=="m";
}

} // namespace

/* parse NMEA GGA ----------------------------------------------------------------
 * returns false for non-GGA sentences, checksum mismatch or malformed fields.
 * fix.valid: quality>=1 and lat/lon present. Empty geoid separation is
 * replaced by the RTKLIB geoid model (embedded EGM96 unless changed with
 * opengeoid()); ellipsoidal height h_ell=h_msl+geoid, ecef from pos2ecef.   */
bool parse_gga(const std::string &sentence, GgaFix &fix)
{
    static const std::string blanks(" \t\r\n\0",5);

    fix=GgaFix();
    size_t a=sentence.find_first_not_of(blanks);
    if (a==std::string::npos) return false;
    size_t b=sentence.find_last_not_of(blanks);
    std::string s=sentence.substr(a,b-a+1);
    if (s.size()<7||s[0]!='$') return false;

    /* checksum: XOR of the characters between '$' and '*' */
    std::string body;
    size_t star=s.find('*');
    if (star!=std::string::npos) {
        if (star+3!=s.size()) return false;
        int h1=hexval(s[star+1]),h2=hexval(s[star+2]);
        if (h1<0||h2<0) return false;
        unsigned cs=0;
        for (size_t i=1;i<star;i++) cs^=(unsigned char)s[i];
        if (cs!=(unsigned)(h1*16+h2)) return false;
        body=s.substr(1,star-1);
    }
    else body=s.substr(1);

    std::vector<std::string> f=split_fields(body);
    if (f.size()<10||f[0].size()!=5||f[0].compare(2,3,"GGA")) return false;
    if (f[0][0]<'A'||f[0][0]>'Z'||f[0][1]<'A'||f[0][1]>'Z') return false;

    /* time, position */
    if (!f[1].empty()&&!to_tod(f[1],fix.tod)) return false;
    bool havelat=!f[2].empty(),havelon=!f[4].empty();
    if (havelat&&!to_angle(f[2],f[3],90.0,'N','S',fix.lat)) return false;
    if (havelon&&!to_angle(f[4],f[5],180.0,'E','W',fix.lon)) return false;

    /* quality, satellites, hdop, altitude above MSL */
    if (!f[6].empty()&&!to_int(f[6],fix.quality)) return false;
    if (!f[7].empty()&&!to_int(f[7],fix.nsat)) return false;
    if (!f[8].empty()&&(!to_num(f[8],fix.hdop)||fix.hdop<0.0)) return false;
    if (!f[9].empty()&&!to_num(f[9],fix.h_msl)) return false;
    if (!is_metre(f,10)||!is_metre(f,12)) return false;

    /* geoid separation: from GGA, else geoid model */
    bool havepos=havelat&&havelon;
    if (f.size()>11&&!f[11].empty()) {
        if (!to_num(f[11],fix.geoid)) return false;
        fix.geoid_from_gga=true;
    }
    else if (havepos) {
        double pos[3]={fix.lat,fix.lon,0.0};
        fix.geoid=geoidh(pos);
    }
    fix.h_ell=fix.h_msl+fix.geoid;
    if (havepos) {
        double pos[3]={fix.lat,fix.lon,fix.h_ell};
        pos2ecef(pos,fix.ecef);
    }
    fix.valid=fix.quality>=1&&havepos;
    return true;
}

/* base64 (RFC 4648) ---------------------------------------------------------------*/
static const char b64_tbl[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

std::string base64_encode(const std::string &in)
{
    std::string out;
    size_t i=0,n=in.size();

    out.reserve((n+2)/3*4);
    for (;i+2<n;i+=3) {
        uint32_t v=(uint32_t)(uint8_t)in[i]<<16|(uint32_t)(uint8_t)in[i+1]<<8|(uint8_t)in[i+2];
        out+=b64_tbl[v>>18&63]; out+=b64_tbl[v>>12&63];
        out+=b64_tbl[v>>6&63];  out+=b64_tbl[v&63];
    }
    if (i<n) {
        uint32_t v=(uint32_t)(uint8_t)in[i]<<16|(i+1<n?(uint32_t)(uint8_t)in[i+1]<<8:0u);
        out+=b64_tbl[v>>18&63]; out+=b64_tbl[v>>12&63];
        out+=i+1<n?b64_tbl[v>>6&63]:'=';
        out+='=';
    }
    return out;
}
/* decode: whitespace ignored, stops at '=', missing padding accepted, URL-safe
 * alphabet accepted. returns empty string on an invalid character          */
std::string base64_decode(const std::string &in)
{
    std::string out;
    uint32_t acc=0;
    int nb=0,v;

    for (char c:in) {
        if      (c>='A'&&c<='Z') v=c-'A';
        else if (c>='a'&&c<='z') v=c-'a'+26;
        else if (c>='0'&&c<='9') v=c-'0'+52;
        else if (c=='+'||c=='-') v=62;
        else if (c=='/'||c=='_') v=63;
        else if (c=='=') break;
        else if (c==' '||c=='\t'||c=='\r'||c=='\n') continue;
        else return std::string();
        acc=(acc<<6|(uint32_t)v)&0xFFFFFF;
        nb+=6;
        if (nb>=8) {nb-=8; out+=(char)(acc>>nb&0xFF);}
    }
    return out;
}

} // namespace nrtk
