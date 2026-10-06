/*------------------------------------------------------------------------------
 * rtcm_out.c : RTCM 3 output (MSM4/MSM7 split to <=64 cells, 1006, 1033,
 *              ephemeris pass-through)
 *----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>
#include "nrtk/rtcm_out.h"

#define RTCM3PREAMB 0xD3

extern int nrtk_rtcmout_init(nrtk_rtcmout_t *o, int staid, int msm,
                             const double *pos, int non_physical)
{
    memset(o,0,sizeof(*o));
    if (!(o->rtcm=(rtcm_t *)calloc(1,sizeof(rtcm_t)))) return 0;
    if (!init_rtcm(o->rtcm)) {free(o->rtcm); o->rtcm=NULL; return 0;}
    o->msm=msm==4?4:7;
    o->staid=staid&0xFFF;
    o->rtcm->staid=o->staid;
    if (pos) memcpy(o->pos,pos,sizeof(o->pos));
    o->non_physical=non_physical;
    o->itrf_year=0;
    snprintf(o->antdes,sizeof(o->antdes),"ADVNULLANTENNA");
    snprintf(o->rectype,sizeof(o->rectype),"NRTK NETWORK RTK");
    snprintf(o->recver,sizeof(o->recver),"%s",NRTK_VERSION);
    o->sysmask=(1<<NRTK_GPS)|(1<<NRTK_GAL)|(1<<NRTK_BDS)|(1<<NRTK_QZS);
    return 1;
}
extern void nrtk_rtcmout_free(nrtk_rtcmout_t *o)
{
    if (o->rtcm) {free_rtcm(o->rtcm); free(o->rtcm);}
    o->rtcm=NULL;
}
extern void nrtk_rtcmout_restart(nrtk_rtcmout_t *o, int staid, const double *pos)
{
    int i,j;
    o->staid=staid&0xFFF;
    o->rtcm->staid=o->staid;
    if (pos) memcpy(o->pos,pos,sizeof(o->pos));
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ+NEXOBS;j++) {
        o->rtcm->lltime[i][j].time=0; o->rtcm->lltime[i][j].sec=0.0;
        o->rtcm->cp[i][j]=0.0;
    }
}
/* bit helpers (RTKLIB setbitu/setbits are exported; 38-bit signed is not) ---*/
static void set38bits(uint8_t *buff, int pos, double value)
{
    int word_h=(int)floor(value/64.0);
    uint32_t word_l=(uint32_t)(value-word_h*64.0);
    setbits(buff,pos,32,word_h);
    setbitu(buff,pos+32,6,word_l);
}
static int frame(uint8_t *msg, int nbit, uint8_t *buf, int size)
{
    int i,len;
    uint32_t crc;
    for (i=nbit;i%8;i++) setbitu(msg,i,1,0);
    len=i/8;
    if (len-3>1023||len+3>size) return 0;
    setbitu(msg,0,8,RTCM3PREAMB);
    setbitu(msg,8,6,0);
    setbitu(msg,14,10,len-3);
    crc=rtk_crc24q(msg,len);
    setbitu(msg,len*8,24,crc);
    memcpy(buf,msg,len+3);
    return len+3;
}
extern int nrtk_rtcm_encode_1006(int staid, const double *pos, double hgt,
                                 int non_physical, int itrf_year, int sysmask,
                                 uint8_t *buf, int size)
{
    uint8_t msg[64]={0};
    int i=24,h=0;

    if (0.0<=hgt&&hgt<=6.5535) h=(int)floor(hgt/0.0001+0.5);
    setbitu(msg,i,12,1006 ); i+=12;              /* DF002 message number */
    setbitu(msg,i,12,staid); i+=12;              /* DF003 station id */
    setbitu(msg,i, 6,itrf_year&0x3F); i+=6;      /* DF021 ITRF realization year */
    setbitu(msg,i, 1,(sysmask&(1<<NRTK_GPS))?1:0); i+=1; /* DF022 GPS */
    setbitu(msg,i, 1,0); i+=1;                   /* DF023 GLONASS */
    setbitu(msg,i, 1,(sysmask&(1<<NRTK_GAL))?1:0); i+=1; /* DF024 Galileo */
    setbitu(msg,i, 1,non_physical?1:0); i+=1;    /* DF141 reference-station indicator */
    set38bits(msg,i,pos[0]/0.0001); i+=38;       /* DF025 ARP ECEF-X */
    setbitu(msg,i, 1,1); i+=1;                   /* DF142 single receiver oscillator */
    setbitu(msg,i, 1,0); i+=1;                   /* DF001 reserved */
    set38bits(msg,i,pos[1]/0.0001); i+=38;       /* DF026 ARP ECEF-Y */
    setbitu(msg,i, 2,0); i+=2;                   /* DF364 quarter cycle indicator */
    set38bits(msg,i,pos[2]/0.0001); i+=38;       /* DF027 ARP ECEF-Z */
    setbitu(msg,i,16,h); i+=16;                  /* DF028 antenna height */
    return frame(msg,i,buf,size);
}
static int append_rtcm(rtcm_t *r, int type, int sync, uint8_t *buf, int pos, int size)
{
    if (!gen_rtcm3(r,type,0,sync)) return pos;
    if (pos+r->nbyte>size) return pos;
    memcpy(buf+pos,r->buff,r->nbyte);
    return pos+r->nbyte;
}
extern int nrtk_rtcmout_station(nrtk_rtcmout_t *o, uint8_t *buf, int size)
{
    int n;
    n=nrtk_rtcm_encode_1006(o->staid,o->pos,o->hgt,o->non_physical,o->itrf_year,
                            o->sysmask,buf,size);
    snprintf(o->rtcm->sta.antdes,sizeof(o->rtcm->sta.antdes),"%s",o->antdes);
    o->rtcm->sta.antsno[0]='\0';
    o->rtcm->sta.antsetup=0;
    snprintf(o->rtcm->sta.rectype,sizeof(o->rtcm->sta.rectype),"%s",o->rectype);
    snprintf(o->rtcm->sta.recver,sizeof(o->rtcm->sta.recver),"%s",o->recver);
    snprintf(o->rtcm->sta.recsno,sizeof(o->rtcm->sta.recsno),"%s",o->recsno);
    o->rtcm->staid=o->staid;
    return append_rtcm(o->rtcm,1033,0,buf,n,size);
}
static int msm_type(int sysidx, int msm)
{
    static const int base[NRTK_NSYS]={1070,1090,1120,1110};
    return base[sysidx]+msm;
}
/* number of distinct MSM signals among obs */
static int nsig_of(const obsd_t *obs, int n)
{
    uint8_t seen[MAXCODE+1]={0};
    int i,j,k=0;
    for (i=0;i<n;i++) for (j=0;j<NFREQ+NEXOBS;j++) {
        uint8_t c=obs[i].code[j];
        if (c==CODE_NONE||c>MAXCODE||seen[c]) continue;
        seen[c]=1; k++;
    }
    return k;
}
extern int nrtk_rtcmout_obs(nrtk_rtcmout_t *o, gtime_t t, const obsd_t *obs, int n,
                            uint8_t *buf, int size)
{
    obsd_t *grp=(obsd_t *)malloc(sizeof(obsd_t)*(n>0?n:1));
    int sys,i,ng,start,pos=0,nmsg=0,imsg=0,pass;
    struct { int sys; int first,count; } msgs[64];

    if (!grp) return 0;
    /* plan messages: per system, satellite groups with nsat*nsig <= 63.
       RTCM allows 64 cells, but RTKLIB-EX 2.5.1 gen_msm_sig() skips the data
       of cell 64 (test "cell>=64" on 1-based cell numbers). */
    for (sys=0;sys<NRTK_NSYS;sys++) {
        if (!(o->sysmask&(1<<sys))) continue;
        for (i=ng=0;i<n;i++) {
            if (nrtk_sysidx(satsys(obs[i].sat,NULL))==sys) grp[ng++]=obs[i];
        }
        for (start=0;start<ng&&nmsg<64;) {
            int cnt=0;
            while (start+cnt<ng) {
                int ns=cnt+1,nsig=nsig_of(grp+start,ns);
                if (ns*nsig>63) break;
                cnt++;
            }
            if (cnt==0) cnt=1;
            msgs[nmsg].sys=sys; msgs[nmsg].first=start; msgs[nmsg].count=cnt; nmsg++;
            start+=cnt;
        }
    }
    /* encode (multiple message bit set except on the last message) */
    for (pass=0;pass<1;pass++) {   /* single pass */
        for (imsg=0;imsg<nmsg;imsg++) {
            int s=msgs[imsg].sys;
            for (i=ng=0;i<n;i++) {
                if (nrtk_sysidx(satsys(obs[i].sat,NULL))==s) grp[ng++]=obs[i];
            }
            o->rtcm->time=t;
            o->rtcm->staid=o->staid;
            o->rtcm->obs.n=msgs[imsg].count;
            memcpy(o->rtcm->obs.data,grp+msgs[imsg].first,sizeof(obsd_t)*msgs[imsg].count);
            pos=append_rtcm(o->rtcm,msm_type(s,o->msm),imsg<nmsg-1,buf,pos,size);
        }
    }
    o->rtcm->obs.n=0;
    free(grp);
    return pos;
}
extern int nrtk_rtcmout_eph(nrtk_rtcmout_t *o, const nav_t *nav, int sat,
                            uint8_t *buf, int size)
{
    int i,best=-1,sys=satsys(sat,NULL),type;

    for (i=0;i<nav->n;i++) {  /* newest ephemeris of the satellite */
        if (nav->eph[i].sat!=sat) continue;
        if (best<0||timediff(nav->eph[i].toe,nav->eph[best].toe)>0.0) best=i;
    }
    if (best<0) return 0;
    switch (sys) {
        case SYS_GPS: type=1019; break;
        case SYS_GAL: type=1046; break;   /* I/NAV */
        case SYS_CMP: type=1042; break;
        case SYS_QZS: type=1044; break;
        default: return 0;
    }
    o->rtcm->ephsat=sat;
    o->rtcm->ephset=0;
    o->rtcm->nav.eph[sat-1]=nav->eph[best];
    return append_rtcm(o->rtcm,type,0,buf,0,size);
}
