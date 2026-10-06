/*------------------------------------------------------------------------------
 * nrtk-sim : synthetic multi-station GNSS observation simulator
 *
 * Generates RINEX 3.04 observation files for a reference station network,
 * a mixed RINEX 3.04 navigation file, a synthetic ANTEX file and an nrtk
 * network configuration, from a truth model:
 *   - Keplerian broadcast ephemerides (GPS, Galileo, BeiDou MEO/IGSO, QZSS QZO)
 *     with a smooth radial/along/cross orbit error for the true orbit
 *   - thin-shell ionosphere: diurnal background, equatorial anomaly crest,
 *     travelling ionospheric disturbance, small-scale structure
 *   - troposphere: Saastamoinen ZHD with height, ZWD with height, horizontal
 *     field, random walk and gradients, mapped with NMF (processing uses GMF)
 *   - receiver antenna PCO/PCV per frequency (written to the ANTEX file),
 *     ARP eccentricity, solid earth tide, receiver clocks (steered or
 *     free-running with millisecond jumps), code/phase biases, multipath,
 *     white noise, cycle slips
 * The truth is self-consistent with RTKLIB-EX conventions (eph2pos, tidedisp)
 * so that leave-one-out validation measures interpolation errors only.
 *----------------------------------------------------------------------------*/
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <yaml-cpp/yaml.h>
#include "nrtk/common.h"
#include "nrtk/antenna.h"

namespace {

const double RE_ION = 6371000.0;
/* constants as in RTKLIB-EX ephemeris.c (not exported by rtklib.h) */
const double MU_GPS=3.9860050E14, MU_GAL=3.986004418E14, MU_CMP=3.986004418E14;
const double OMGE_GAL=7.2921151467E-5, OMGE_CMP=7.292115E-5;
inline double SQR(double x) {return x*x;}
const double K_ION  = 40.3E16/(FREQL1*FREQL1); /* m of L1 delay per TECU */

struct SigDef { const char *code; int sysmask_prn_lo, prn_hi; };

struct Station {
    std::string id, antenna;
    double llh[3]={0}, pos[3]={0}, antdel[3]={0};
    bool free_clock=false;
    double mp_scale=1.0;
    double slip_scale=1.0;      /* multiplies the slip rates for this station (0: none) */
    double zwd_rw=0.0, pres_off=0.0;
    double clk=0.0, clk_drift=0.0;
    double bias_code[NRTK_NSYS][8]={{0}}, bias_phase[NRTK_NSYS][8]={{0}};
    const nrtk_ant_t *ant=nullptr;
};

struct SatTrack {               /* per station per satellite state */
    bool  active=false;
    double amb[NRTK_MAXSIG]={0};/* integer ambiguities (cycles) */
    double mp_ph[2]={0}, mp_per[2]={0};
};

struct Scenario {
    gtime_t t0;
    double dur=3600, tint=1.0;
    unsigned seed=1;
    int sysmask=0xF;
    std::string outdir="data/sim/default";
    std::vector<Station> sta;
    /* ionosphere */
    double vtec_day=50, vtec_night=10, crest_lat=23.0, crest_amp=20, crest_w=3.0;
    double tid_amp=0.3, tid_wl=300, tid_per=40, tid_az=160, shell_h=350;
    double ss_amp=0.05, storm=0.0;
    /* troposphere */
    double zwd_sl=0.28, zwd_hs=2000, zwd_gn=0.01, zwd_ge=0.005, zwd_rw=0.01;
    double zwd_rw_local=0.004, grad=0.0005;
    /* orbit error amplitudes RAC (m) */
    double orb[3]={0.3,1.0,0.6};
    double noise_ph=0.0015, noise_code=0.15, slip_rate=2E-5, elmin=5.0;
    double slip_unflagged_rate=0.0; /* small slips without LLI (receiver misses them) */
    int    slip_unflagged_max=2;    /* max |cycles| of an unflagged slip */
    double slip_common_rate=0.0;    /* unflagged slips of the same size on all signals of
                                       a satellite (invisible to geometry-free and MW) */
    double ref_lat=24.8, ref_lon=121.1;
    double truth_int=10.0;      /* truth output interval (s), 0: off */
};

std::mt19937_64 rng;
double gauss() { static std::normal_distribution<double> n(0.0,1.0); return n(rng); }
double unif()  { static std::uniform_real_distribution<double> u(0.0,1.0); return u(rng); }

/* ---------------------------------------------------------------- orbits ---*/
struct OrbitDef { int sat; double A,e,i,O0,w,M0; double phase[3]; };

void kepler_set(std::vector<OrbitDef> &orb, int sys, int prn, double A,
                double e, double inc, double O, double w, double M)
{
    OrbitDef d;
    d.sat=satno(sys,prn); if (!d.sat) return;
    d.A=A; d.e=e; d.i=inc*D2R; d.O0=O*D2R; d.w=w*D2R; d.M0=M*D2R;
    for (int k=0;k<3;k++) d.phase[k]=unif()*2*PI;
    orb.push_back(d);
}
std::vector<OrbitDef> constellation(int sysmask)
{
    std::vector<OrbitDef> o;
    if (sysmask&(1<<NRTK_GPS)) {      /* 6 planes x 5 */
        for (int p=0;p<6;p++) for (int s=0;s<5;s++)
            kepler_set(o,SYS_GPS,p*5+s+1,26559710.0,0.005+0.01*unif(),55.0,
                       p*60.0+7.0,unif()*360.0,s*72.0+p*12.0);
    }
    if (sysmask&(1<<NRTK_GAL)) {      /* Walker 24/3/1 */
        for (int p=0;p<3;p++) for (int s=0;s<8;s++)
            kepler_set(o,SYS_GAL,p*8+s+1,29599801.0,0.0005,56.0,p*120.0+25.0,0.0,
                       s*45.0+p*15.0);
    }
    if (sysmask&(1<<NRTK_BDS)) {      /* MEO Walker 24/3/1 + 3 IGSO */
        for (int p=0;p<3;p++) for (int s=0;s<8;s++)
            kepler_set(o,SYS_CMP,19+p*8+s,27906100.0,0.001,55.0,p*120.0+50.0,0.0,
                       s*45.0+p*15.0);
        for (int s=0;s<3;s++)         /* node longitude 118E (O0 set later) */
            kepler_set(o,SYS_CMP,6+s,42164000.0,0.002,55.0,118.0,180.0,s*120.0);
    }
    if (sysmask&(1<<NRTK_QZS)) {      /* QZO, figure-8 around 135E */
        int prns[3]={193,194,196};
        for (int s=0;s<3;s++)
            kepler_set(o,SYS_QZS,prns[s],42164000.0,0.075,41.0,135.0,270.0,s*120.0);
    }
    return o;
}
double sys_mu(int sys)   { return sys==SYS_GAL?MU_GAL:(sys==SYS_CMP?MU_CMP:MU_GPS); }
double sys_omge(int sys) { return sys==SYS_GAL?OMGE_GAL:(sys==SYS_CMP?OMGE_CMP:OMGE); }

/* build broadcast ephemeris for orbit d at toe ------------------------------*/
eph_t make_eph(const OrbitDef &d, gtime_t toe, gtime_t tref, int iode, double af0,
               double af1)
{
    int sys=satsys(d.sat,NULL), week;
    double n=sqrt(sys_mu(sys)/(d.A*d.A*d.A)), omge=sys_omge(sys);
    eph_t e={0};
    e.sat=d.sat; e.iode=iode; e.iodc=iode; e.sva=0; e.svh=0;
    e.A=d.A; e.e=d.e; e.i0=d.i; e.omg=d.w; e.deln=0.0; e.idot=0.0;
    e.OMGd=(sys==SYS_GPS||sys==SYS_GAL||sys==SYS_CMP)&&d.A<4E7?-2.6E-9:0.0;
    double tk=timediff(toe,tref);
    e.M0=d.M0+n*tk;
    /* BDS: toe and week in BDT (week counted from 2006), others GPS week */
    double toes=sys==SYS_CMP?time2bdt(gpst2bdt(toe),&week):time2gpst(toe,&week);
    /* node: O(t) = OMG0 + (OMGd-OMGE)*(t-toe) - OMGE*toes  must equal
       O_ref(t) = O0 + OMGd*(t-tref) - OMGE*(t-tref) - OMGE*toes_ref         */
    double O0=d.O0;
    if (d.A>4E7) { /* geosynchronous: node longitude O0 at the crossing */
        O0=d.O0-d.M0-d.w;
    }
    e.OMG0=O0+e.OMGd*tk-omge*tk+omge*toes; /* tref: reference with toes_ref=0 */
    e.OMG0=fmod(e.OMG0,2*PI);
    e.toes=toes; e.week=week;
    e.toe=toe; e.toc=toe; e.ttr=timeadd(toe,-3600.0);
    e.f0=af0; e.f1=af1; e.f2=0.0; e.fit=4.0;
    if (sys==SYS_GAL) e.code=(1<<0)|(1<<2)|(1<<9); /* I/NAV */
    if (sys==SYS_CMP) e.flag=1;
    return e;
}

/* true satellite position: broadcast + RAC error ----------------------------*/
bool truth_satpos(gtime_t t, int sat, const nav_t *nav,
                  const std::map<int,OrbitDef> &orb, const double *amp,
                  double *rs, double *dts)
{
    double var; int svh;
    if (!satpos(t,t,sat,EPHOPT_BRDC,nav,rs,dts,&var,&svh)) return false;
    auto it=orb.find(sat);
    if (it==orb.end()) return true;
    const OrbitDef &d=it->second;
    double r=norm(rs,3),er[3],ec[3],ea[3],h[3];
    for (int i=0;i<3;i++) er[i]=rs[i]/r;
    cross3(rs,rs+3,h);
    double hn=norm(h,3);
    for (int i=0;i<3;i++) ec[i]=h[i]/hn;
    cross3(ec,er,ea);
    double T=2*PI*sqrt(d.A*d.A*d.A/sys_mu(satsys(sat,NULL)));
    double tt=time2gpst(t,NULL);
    double dr=amp[0]*sin(2*PI*tt/T+d.phase[0]);
    double da=amp[1]*sin(2*PI*tt/T+d.phase[1])+0.3*amp[1]*sin(2*PI*tt/43200.0+d.phase[2]);
    double dc=amp[2]*sin(2*PI*tt/T+d.phase[2]);
    for (int i=0;i<3;i++) rs[i]+=dr*er[i]+da*ea[i]+dc*ec[i];
    return true;
}

/* ------------------------------------------------------------ ionosphere ---*/
double vtec(const Scenario &sc, double lat, double lon, gtime_t t)
{
    double ep[6]; time2epoch(gpst2utc(t),ep);
    double ut=ep[3]+ep[4]/60.0+ep[5]/3600.0;
    double lt=fmod(ut+lon*R2D/15.0+48.0,24.0);
    double day=cos(2*PI*(lt-14.0)/24.0); day=day>0?pow(day,1.5):0.0;
    double crest_t=exp(-pow((lt-16.0)/4.0,2)); /* afternoon / post-sunset crest */
    double latd=lat*R2D, lond=lon*R2D;
    double v=sc.vtec_night+(sc.vtec_day-sc.vtec_night)*day;
    v+=sc.crest_amp*crest_t*exp(-pow((latd-sc.crest_lat)/sc.crest_w,2));
    /* TID: plane wave */
    double x=(lond-sc.ref_lon)*111.32*cos(sc.ref_lat*D2R), y=(latd-sc.ref_lat)*110.57;
    double d=x*sin(sc.tid_az*D2R)+y*cos(sc.tid_az*D2R);
    double tt=time2gpst(t,NULL);
    v+=sc.tid_amp*sin(2*PI*(d/sc.tid_wl-tt/(sc.tid_per*60.0)));
    /* small scale structure */
    v+=sc.ss_amp*(sin(2*PI*(x/70.0+tt/900.0))*cos(2*PI*(y/55.0-tt/1300.0))+
                  0.7*sin(2*PI*(x/40.0-y/45.0+tt/700.0)));
    /* disturbed: travelling depletion band after sunset */
    if (sc.storm>0.0) {
        double band=exp(-pow((x-0.05*(tt-3600*((int)(tt/3600))))/(25.0),2));
        v-=sc.storm*band*exp(-pow((lt-21.0)/2.0,2));
    }
    return v<0.5?0.5:v;
}
double iono_l1(const Scenario &sc, const double *pos, const double *azel, gtime_t t)
{
    double H=sc.shell_h*1E3, el=azel[1], az=azel[0];
    double rp=RE_ION/(RE_ION+H)*cos(el);
    double psi=PI/2-el-asin(rp);
    double lat=asin(sin(pos[0])*cos(psi)+cos(pos[0])*sin(psi)*cos(az));
    double lon=pos[1]+asin(sin(psi)*sin(az)/cos(lat));
    double F=1.0/sqrt(1.0-rp*rp);
    return K_ION*F*vtec(sc,lat,lon,t);
}

/* ---------------------------------------------------------- troposphere ---*/
double trop_true(const Scenario &sc, Station &s, const double *pos,
                 const double *azel, gtime_t t)
{
    double hgt=pos[2]-geoidh(pos); if (hgt<0) hgt=0;
    double pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568)+s.pres_off;
    double zhd=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3);
    double x=(pos[1]*R2D-sc.ref_lon)*111.32*cos(sc.ref_lat*D2R)/100.0;
    double y=(pos[0]*R2D-sc.ref_lat)*110.57/100.0;
    double zwd=sc.zwd_sl*exp(-hgt/sc.zwd_hs)+sc.zwd_ge*x+sc.zwd_gn*y+s.zwd_rw;
    double tt=time2gpst(t,NULL);
    zwd+=0.01*sin(2*PI*tt/14400.0+0.3*x);         /* common slow variation */
    if (zwd<0.01) zwd=0.01;
    double mw,mh=tropmapf(t,pos,azel,&mw);        /* NMF */
    double cotel=1.0/tan(azel[1]);
    double gn=sc.grad*cos(2*PI*tt/20000.0), ge=sc.grad*sin(2*PI*tt/26000.0);
    return zhd*mh+zwd*mw+mw*cotel*(gn*cos(azel[0])+ge*sin(azel[0]));
}

/* ----------------------------------------------------------- signals -------*/
struct SysSig { int sys; std::vector<std::string> codes; };

std::vector<std::string> sat_signals(int sat)
{
    int prn, sys=satsys(sat,&prn);
    std::vector<std::string> v;
    switch (sys) {
        case SYS_GPS:
            v={"1C","2W"};
            if (prn<=20) v.push_back("2L");
            if (prn<=14) v.push_back("5Q");
            break;
        case SYS_GAL: v={"1C","5Q","7Q"}; break;
        case SYS_CMP:
            v={"2I","6I"};
            if (prn>=19&&prn<=46) {v.push_back("1P"); v.push_back("5P");}
            break;
        case SYS_QZS: v={"1C","2L","5Q"}; break;
    }
    return v;
}

/* ------------------------------------------------------------- ANTEX -------*/
void write_antex(const std::string &file, const std::vector<std::string> &types)
{
    FILE *fp=fopen(file.c_str(),"w");
    if (!fp) {fprintf(stderr,"cannot write %s\n",file.c_str()); exit(1);}
    fprintf(fp,"%8.1f%12s%-1s%-19s%-20s%-20s\n",1.4,"","M","","","ANTEX VERSION / SYST");
    fprintf(fp,"%-1s%-59s%-20s\n","A","","PCV TYPE / REFANT");
    fprintf(fp,"%-60s%-20s\n","synthetic antennas generated by nrtk-sim","COMMENT");
    fprintf(fp,"%-60s%-20s\n","","END OF HEADER");
    const char *fr[]={"G01","G02","G05","E01","E05","E07","E06","E08","C01","C02",
                      "C05","C06","C07","J01","J02","J05"};
    for (size_t a=0;a<types.size();a++) {
        double sc=1.0+0.25*a;
        fprintf(fp,"%60s%-20s\n","","START OF ANTENNA");
        char tr[64]; snprintf(tr,sizeof(tr),"%-20s",types[a].c_str());
        fprintf(fp,"%-20.20s%-20s%-20s%-20s\n",tr,"","","TYPE / SERIAL NO");
        fprintf(fp,"%-20s%-20s%10d%10s%-20s\n","ROBOT","nrtk-sim",0,"","METH / BY / # / DATE");
        fprintf(fp,"%8.1f%52s%-20s\n",5.0,"","DAZI");
        fprintf(fp,"%8.1f%6.1f%6.1f%40s%-20s\n",0.0,90.0,5.0,"","ZEN1 / ZEN2 / DZEN");
        fprintf(fp,"%6d%54s%-20s\n",(int)(sizeof(fr)/sizeof(fr[0])),"","# OF FREQUENCIES");
        for (size_t f=0;f<sizeof(fr)/sizeof(fr[0]);f++) {
            double band=(fr[f][2]=='1'||(fr[f][2]=='2'&&fr[f][0]=='C'))?1.0:
                        (fr[f][2]=='2'?0.85:0.75);
            double up=(62.0+10.0*band+3.0*f/16.0)*sc, n0=(0.8+0.3*f/16.0)*sc, e0=(-0.4+0.2*band)*sc;
            fprintf(fp,"   %-3s%54s%-20s\n",fr[f],"","START OF FREQUENCY");
            fprintf(fp,"%10.2f%10.2f%10.2f%30s%-20s\n",n0,e0,up,"","NORTH / EAST / UP");
            fprintf(fp,"   NOAZI");
            double pcv[19];
            for (int z=0;z<19;z++) {
                double zen=z*5.0*D2R;
                pcv[z]=sc*band*(-1.5*sin(2*zen)+4.0*(1-cos(zen))-2.0*pow(sin(zen),4));
                fprintf(fp,"%8.2f",pcv[z]);
            }
            fprintf(fp,"\n");
            for (int ia=0;ia<=72;ia++) {
                double az=ia*5.0*D2R;
                fprintf(fp,"%8.1f",ia*5.0);
                for (int z=0;z<19;z++) {
                    double zen=z*5.0*D2R;
                    fprintf(fp,"%8.2f",pcv[z]+0.8*sc*sin(zen)*cos(az-0.6*f)+
                            0.4*sc*sin(2*zen)*sin(2*az));
                }
                fprintf(fp,"\n");
            }
            fprintf(fp,"   %-3s%54s%-20s\n",fr[f],"","END OF FREQUENCY");
        }
        fprintf(fp,"%60s%-20s\n","","END OF ANTENNA");
    }
    fclose(fp);
}

/* --------------------------------------------------------- config parse ----*/
gtime_t parse_time(const std::string &s)
{
    double ep[6]={0};
    if (sscanf(s.c_str(),"%lf-%lf-%lf %lf:%lf:%lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<3 &&
        sscanf(s.c_str(),"%lf/%lf/%lf %lf:%lf:%lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<3) {
        fprintf(stderr,"bad time: %s\n",s.c_str()); exit(1);
    }
    return epoch2time(ep);
}
template<class T> T get(const YAML::Node &n, const char *k, T def)
{
    return n&&n[k]?n[k].as<T>():def;
}
Scenario load(const std::string &file)
{
    Scenario sc;
    YAML::Node y=YAML::LoadFile(file);
    sc.t0=parse_time(get<std::string>(y,"start","2026-09-20 04:00:00"));
    sc.dur=get<double>(y,"duration_s",3600);
    sc.tint=get<double>(y,"interval_s",1.0);
    sc.seed=get<unsigned>(y,"seed",1);
    sc.outdir=get<std::string>(y,"output_dir",sc.outdir);
    sc.elmin=get<double>(y,"elevation_min_deg",5.0);
    if (y["systems"]) {
        std::string s;
        for (auto v:y["systems"]) s+=v.as<std::string>()+",";
        sc.sysmask=nrtk_sysmask(s.c_str());
    }
    YAML::Node io=y["iono"];
    sc.vtec_day=get(io,"vtec_day",sc.vtec_day); sc.vtec_night=get(io,"vtec_night",sc.vtec_night);
    sc.crest_lat=get(io,"crest_lat",sc.crest_lat); sc.crest_amp=get(io,"crest_amp",sc.crest_amp);
    sc.crest_w=get(io,"crest_width_deg",sc.crest_w); sc.tid_amp=get(io,"tid_amp",sc.tid_amp);
    sc.tid_wl=get(io,"tid_wavelength_km",sc.tid_wl); sc.tid_per=get(io,"tid_period_min",sc.tid_per);
    sc.tid_az=get(io,"tid_azimuth_deg",sc.tid_az); sc.shell_h=get(io,"shell_height_km",sc.shell_h);
    sc.ss_amp=get(io,"small_scale_amp",sc.ss_amp); sc.storm=get(io,"depletion_amp",sc.storm);
    YAML::Node tr=y["trop"];
    sc.zwd_sl=get(tr,"zwd_sealevel",sc.zwd_sl); sc.zwd_hs=get(tr,"zwd_scale_height_m",sc.zwd_hs);
    sc.zwd_gn=get(tr,"zwd_north_m_per_100km",sc.zwd_gn); sc.zwd_ge=get(tr,"zwd_east_m_per_100km",sc.zwd_ge);
    sc.zwd_rw=get(tr,"zwd_rw_m_per_sqrt_h",sc.zwd_rw); sc.zwd_rw_local=get(tr,"zwd_local_rw_m_per_sqrt_h",sc.zwd_rw_local);
    sc.grad=get(tr,"gradient_m",sc.grad);
    if (y["orbit_error_rac_m"]) for (int i=0;i<3;i++) sc.orb[i]=y["orbit_error_rac_m"][i].as<double>();
    YAML::Node no=y["noise"];
    sc.noise_ph=get(no,"phase_m",sc.noise_ph); sc.noise_code=get(no,"code_m",sc.noise_code);
    sc.slip_rate=get(y,"slip_rate",sc.slip_rate);
    sc.slip_unflagged_rate=get(y,"slip_unflagged_rate",sc.slip_unflagged_rate);
    sc.slip_unflagged_max=get(y,"slip_unflagged_max",sc.slip_unflagged_max);
    sc.slip_common_rate=get(y,"slip_common_rate",sc.slip_common_rate);
    sc.ref_lat=get(y,"ref_lat",sc.ref_lat); sc.ref_lon=get(y,"ref_lon",sc.ref_lon);
    sc.truth_int=get(y,"truth_interval_s",sc.truth_int);
    for (auto n:y["stations"]) {
        Station s;
        s.id=n["id"].as<std::string>();
        for (int i=0;i<3;i++) s.llh[i]=n["llh"][i].as<double>();
        s.antenna=get<std::string>(n,"antenna","NRTKSIM_ANT1 NONE");
        s.antdel[2]=get<double>(n,"arp_height_m",0.0);
        s.free_clock=get<std::string>(n,"clock","steered")=="free";
        s.mp_scale=get<double>(n,"multipath",1.0);
        s.slip_scale=get<double>(n,"slip_scale",1.0);
        double pos[3]={s.llh[0]*D2R,s.llh[1]*D2R,s.llh[2]};
        pos2ecef(pos,s.pos);
        sc.sta.push_back(s);
    }
    return sc;
}

void setup_rnxopt(rnxopt_t &opt, const Scenario &sc, const Station &s)
{
    memset(&opt,0,sizeof(opt));
    opt.rnxver=304;
    opt.navsys=SYS_GPS|SYS_GAL|SYS_CMP|SYS_QZS;
    opt.obstype=OBSTYPE_PR|OBSTYPE_CP|OBSTYPE_DOP|OBSTYPE_SNR;
    opt.freqtype=FREQTYPE_ALL;
    opt.tint=sc.tint;
    for (int m=0;m<RNX_NUMSYS;m++) memset(opt.mask[m],'1',MAXCODE);
    strcpy(opt.prog,"nrtk-sim "); strcpy(opt.runby,"nrtk");
    snprintf(opt.marker,sizeof(opt.marker),"%s",s.id.c_str());
    strcpy(opt.markertype,"GEODETIC");
    strcpy(opt.rec[1],"NRTK SIM RECEIVER"); strcpy(opt.rec[2],"1.0");
    char ant[64]; snprintf(ant,sizeof(ant),"%s",s.antenna.c_str());
    char type[32]="",radome[16]="NONE";
    sscanf(ant,"%31s %15s",type,radome);
    snprintf(opt.ant[1],sizeof(opt.ant[1]),"%-16s%-4s",type,radome);
    for (int i=0;i<3;i++) opt.apppos[i]=s.pos[i];
    opt.antdel[0]=s.antdel[2];
    strcpy(opt.comment[0],"synthetic data: nrtk-sim truth model");
    struct {int rs; int sys;} m[]={{RNX_SYS_GPS,SYS_GPS},{RNX_SYS_GAL,SYS_GAL},
                                   {RNX_SYS_CMP,SYS_CMP},{RNX_SYS_QZS,SYS_QZS}};
    const char *sigs[4][6]={{"1C","2W","2L","5Q",0,0},{"1C","5Q","7Q",0,0,0},
                            {"2I","6I","1P","5P",0,0},{"1C","2L","5Q",0,0,0}};
    for (int k=0;k<4;k++) {
        int n=0;
        for (int j=0;sigs[k][j];j++) for (const char *t: {"C","L","D","S"}) {
            snprintf(opt.tobs[m[k].rs][n++],4,"%s%s",t,sigs[k][j]);
        }
        opt.nobs[m[k].rs]=n;
    }
    opt.tstart=sc.t0;
    opt.tend=timeadd(sc.t0,sc.dur-sc.tint);
}

} // namespace

int main(int argc, char **argv)
{
    if (argc<2) {
        fprintf(stderr,"usage: nrtk-sim scenario.yaml\n");
        return 1;
    }
    Scenario sc=load(argv[1]);
    rng.seed(sc.seed);
    for (size_t i=1;i<=sc.outdir.size();i++) { /* mkdir -p */
        if (i==sc.outdir.size()||sc.outdir[i]=='/') mkdir(sc.outdir.substr(0,i).c_str(),0755);
    }

    /* antennas */
    std::vector<std::string> anttypes;
    for (auto &s:sc.sta) {
        bool f=false;
        for (auto &t:anttypes) if (t==s.antenna) f=true;
        if (!f) anttypes.push_back(s.antenna);
    }
    std::string atxfile=sc.outdir+"/nrtksim.atx";
    {
        std::vector<std::string> fmt;
        for (auto &t:anttypes) {
            char type[32]="",radome[16]="NONE",b[32];
            sscanf(t.c_str(),"%31s %15s",type,radome);
            snprintf(b,sizeof(b),"%-16s%-4s",type,radome);
            fmt.push_back(b);
        }
        write_antex(atxfile,fmt);
    }
    nrtk_atx_t atx={0};
    nrtk_atx_read(atxfile.c_str(),&atx);
    for (auto &s:sc.sta) {
        s.ant=nrtk_atx_find(&atx,s.antenna.c_str());
        if (!s.ant) {fprintf(stderr,"antenna %s not found\n",s.antenna.c_str()); return 1;}
    }

    /* constellation and broadcast ephemerides */
    std::vector<OrbitDef> orbs=constellation(sc.sysmask);
    std::map<int,OrbitDef> orbmap;
    nav_t nav={0};
    gtime_t tref=gpst2time(0,0.0);
    {
        double ep[6]; time2epoch(sc.t0,ep);
        int week; time2gpst(sc.t0,&week);
        tref=gpst2time(week,0.0);
    }
    gtime_t ts=timeadd(sc.t0,-7200.0), te=timeadd(sc.t0,sc.dur+7200.0);
    std::map<int,std::pair<double,double>> clk;
    for (auto &d:orbs) {
        orbmap[d.sat]=d;
        clk[d.sat]={(unif()-0.5)*1E-3,(unif()-0.5)*1E-11};
    }
    int iode=1;
    for (auto &d:orbs) {
        int sys=satsys(d.sat,NULL);
        double step=sys==SYS_GAL?3600.0:(sys==SYS_CMP?3600.0:7200.0);
        double t0s=floor(time2gpst(ts,NULL)/step)*step;
        int week; time2gpst(ts,&week);
        for (gtime_t toe=gpst2time(week,t0s);timediff(toe,te)<=0;toe=timeadd(toe,step)) {
            double dtk=timediff(toe,tref);
            eph_t e=make_eph(d,toe,tref,(iode++)%240+1,clk[d.sat].first+clk[d.sat].second*dtk,
                             clk[d.sat].second);
            if (nav.n>=nav.nmax) {
                nav.nmax=nav.nmax?nav.nmax*2:1024;
                nav.eph=(eph_t *)realloc(nav.eph,sizeof(eph_t)*nav.nmax);
            }
            nav.eph[nav.n++]=e;
        }
    }
    /* navigation file */
    {
        rnxopt_t opt;
        memset(&opt,0,sizeof(opt));
        opt.rnxver=304; opt.navsys=SYS_GPS|SYS_GAL|SYS_CMP|SYS_QZS;
        strcpy(opt.prog,"nrtk-sim"); strcpy(opt.runby,"nrtk");
        std::string f=sc.outdir+"/brdc.nav";
        FILE *fp=fopen(f.c_str(),"w");
        outrnxnavh(fp,&opt,&nav);
        for (int i=0;i<nav.n;i++) outrnxnavb(fp,&opt,nav.eph+i);
        fclose(fp);
    }

    /* station setup */
    for (auto &s:sc.sta) {
        for (int k=0;k<NRTK_NSYS;k++) for (int j=0;j<8;j++) {
            s.bias_code[k][j]=gauss()*0.5;
            s.bias_phase[k][j]=unif();
        }
        s.clk=s.free_clock?(unif()-0.5)*1E-3:gauss()*20E-9;
        s.clk_drift=s.free_clock?(unif()+0.5)*3E-8:0.0;
        s.pres_off=gauss()*2.0;
    }
    std::map<int,double[8]> satbias_code, satbias_phase;
    for (auto &d:orbs) for (int j=0;j<8;j++) {
        satbias_code[d.sat][j]=gauss()*1.0; satbias_phase[d.sat][j]=unif();
    }

    /* observation files */
    std::vector<FILE *> fps;
    std::vector<rnxopt_t> opts(sc.sta.size());
    for (size_t k=0;k<sc.sta.size();k++) {
        setup_rnxopt(opts[k],sc,sc.sta[k]);
        std::string f=sc.outdir+"/"+sc.sta[k].id+".rnx";
        FILE *fp=fopen(f.c_str(),"w");
        if (!fp) {fprintf(stderr,"cannot write %s\n",f.c_str()); return 1;}
        outrnxobsh(fp,&opts[k],&nav);
        fps.push_back(fp);
    }
    std::vector<FILE *> ftruth(sc.sta.size(),nullptr);
    if (sc.truth_int>0.0) {
        for (size_t k=0;k<sc.sta.size();k++) {
            ftruth[k]=fopen((sc.outdir+"/truth_"+sc.sta[k].id+".csv").c_str(),"w");
            fprintf(ftruth[k],"tow,sat,az,el,iono_l1,trop,orbit\n");
        }
    }
    std::vector<std::map<int,SatTrack>> track(sc.sta.size());
    std::vector<obsd_t> obs(MAXOBS);
    double zwd_common=0.0;
    long nep=(long)(sc.dur/sc.tint+0.5);
    long nslip=0,nunflagged=0;

    for (long ie=0;ie<nep;ie++) {
        gtime_t t=timeadd(sc.t0,ie*sc.tint);
        zwd_common+=sc.zwd_rw*sqrt(sc.tint/3600.0)*gauss();
        for (size_t k=0;k<sc.sta.size();k++) {
            Station &s=sc.sta[k];
            /* receiver clock */
            if (ie>0) {
                s.clk+=s.clk_drift*sc.tint+(s.free_clock?2E-11:1E-10)*gauss();
                if (!s.free_clock) s.clk*=0.995;
                if (s.free_clock&&fabs(s.clk)>0.5E-3) s.clk-=(s.clk>0?1E-3:-1E-3);
            }
            s.zwd_rw=zwd_common+0.0;
            static std::vector<double> loc(64,0.0);
            loc[k]+=sc.zwd_rw_local*sqrt(sc.tint/3600.0)*gauss();
            s.zwd_rw+=loc[k];

            double r0[3],rr[3],pos[3];
            nrtk_tide_pos(t,s.pos,rr);
            for (int i=0;i<3;i++) r0[i]=s.pos[i];
            ecef2pos(rr,pos);
            gtime_t trx=timeadd(t,-s.clk);
            int n=0;
            for (auto &d:orbs) {
                int sat=d.sat, sys=satsys(sat,NULL), idx=nrtk_sysidx(sys);
                double rs[6],dts[2],tau=0.075,rho=0,dr[3],e[3],azel[2];
                bool ok=true;
                for (int it=0;it<5;it++) {
                    gtime_t ttx=timeadd(trx,-tau);
                    if (!truth_satpos(ttx,sat,&nav,orbmap,sc.orb,rs,dts)) {ok=false; break;}
                    for (int i=0;i<3;i++) dr[i]=rs[i]-rr[i];
                    rho=norm(dr,3);
                    tau=rho/CLIGHT;
                }
                if (!ok) continue;
                for (int i=0;i<3;i++) e[i]=dr[i]/rho;
                double r=rho+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
                double rate=dot3(rs+3,e);
                satazel(pos,e,azel);
                SatTrack &tk=track[k][sat];
                if (azel[1]<sc.elmin*D2R) {tk.active=false; continue;}
                bool rising=!tk.active;
                if (rising) {
                    tk.active=true;
                    for (int j=0;j<2;j++) {
                        tk.mp_ph[j]=unif()*2*PI;
                        tk.mp_per[j]=300.0+600.0*unif();
                    }
                }
                double T=trop_true(sc,s,pos,azel,t);
                double I=iono_l1(sc,pos,azel,t);
                double tt=time2gpst(t,NULL);
                double sinel=sin(azel[1]);
                if (ftruth[k]&&fmod(tt+1E-6,sc.truth_int)<1E-3) {
                    /* broadcast-orbit range at the same transmission time */
                    double rsb[6],dtsb[2],var; int svh;
                    gtime_t ttx=timeadd(trx,-tau);
                    satpos(ttx,ttx,sat,EPHOPT_BRDC,&nav,rsb,dtsb,&var,&svh);
                    double db[3]; for (int q=0;q<3;q++) db[q]=rsb[q]-rr[q];
                    double rb=norm(db,3)+OMGE*(rsb[0]*rr[1]-rsb[1]*rr[0])/CLIGHT;
                    char id[8]; satno2id(sat,id);
                    fprintf(ftruth[k],"%.1f,%s,%.4f,%.4f,%.5f,%.5f,%.5f\n",tt,id,azel[0]*R2D,
                            azel[1]*R2D,I,T,r-rb);
                }
                double mpw=s.mp_scale*exp(-azel[1]*R2D/25.0);
                double mp_l=0.004*mpw*sin(2*PI*tt/tk.mp_per[0]+tk.mp_ph[0]);
                double mp_p=0.35*mpw*sin(2*PI*tt/tk.mp_per[1]+tk.mp_ph[1]);
                std::vector<std::string> sigs=sat_signals(sat);
                obsd_t &o=obs[n];
                memset(&o,0,sizeof(o));
                o.time=t; o.sat=sat; o.rcv=1;
                double common_slip=0.0;
                if (!rising&&unif()<sc.slip_common_rate*s.slip_scale) {
                    double k=floor(unif()*sc.slip_unflagged_max)+1.0;
                    common_slip=unif()<0.5?k:-k; nslip++; nunflagged++;
                }
                for (size_t j=0;j<sigs.size()&&j<NRTK_MAXSIG;j++) {
                    uint8_t code=obs2code(sigs[j].c_str());
                    double f=nrtk_code2freq(sys,code), lam=CLIGHT/f, mu=SQR(FREQL1/f);
                    double dant=nrtk_ant_dr(s.ant,sys,f,azel,s.antdel,1);
                    double geo=r+CLIGHT*(s.clk-dts[0])+T+dant;
                    double bc=s.bias_code[idx][j]-satbias_code[sat][j];
                    double bp=s.bias_phase[idx][j]-satbias_phase[sat][j];
                    double ph_noise=sc.noise_ph*(0.5+0.5/sinel)*gauss();
                    double cd_noise=sc.noise_code*(0.5+0.5/sinel)*gauss();
                    double P=geo+mu*I+bc+mp_p*(j==0?1.0:0.8)+cd_noise;
                    double Lm=geo-mu*I+mp_l*(j==0?1.0:1.2)+ph_noise;
                    if (rising) tk.amb[j]=floor((P-Lm)/lam+0.5)+floor(unif()*20-10);
                    o.LLI[j]=0;
                    if (!rising&&unif()<sc.slip_rate*s.slip_scale) {
                        tk.amb[j]+=floor(unif()*20.0-10.0)+1.0; o.LLI[j]|=LLI_SLIP; nslip++;
                    }
                    else if (!rising&&unif()<sc.slip_unflagged_rate*s.slip_scale) {
                        /* small slip the receiver does not flag: +-1..max cycles */
                        double k=floor(unif()*sc.slip_unflagged_max)+1.0;
                        tk.amb[j]+=unif()<0.5?k:-k; nslip++; nunflagged++;
                    }
                    tk.amb[j]+=common_slip;
                    o.code[j]=code;
                    o.P[j]=P;
                    o.L[j]=Lm/lam+tk.amb[j]+bp;
                    o.D[j]=(float)(-(rate+CLIGHT*(s.clk_drift-dts[1]))/lam);
                    o.SNR[j]=(float)(32.0+16.0*sinel-(j>0?2.0:0.0)+gauss()*0.5);
                }
                n++;
                if (n>=MAXOBS) break;
            }
            if (n>0) {
                std::sort(obs.begin(),obs.begin()+n,
                          [](const obsd_t &a,const obsd_t &b){return a.sat<b.sat;});
                outrnxobsb(fps[k],&opts[k],obs.data(),n,0);
            }
        }
        if (ie%600==0) fprintf(stderr,"\r  epoch %ld/%ld",ie,nep);
    }
    fprintf(stderr,"\r  %ld epochs, %ld slips (%ld unflagged)          \n",nep,nslip,nunflagged);
    for (auto fp:fps) fclose(fp);
    for (auto fp:ftruth) if (fp) fclose(fp);

    /* truth and processing configuration */
    std::string cf=sc.outdir+"/network.yaml";
    FILE *fp=fopen(cf.c_str(),"w");
    double ep[6]; time2epoch(sc.t0,ep);
    fprintf(fp,"# generated by nrtk-sim from %s\n",argv[1]);
    fprintf(fp,"network:\n  name: SIM\n  frame: ITRF2020\n  epoch: %.4f\n  max_baseline_km: 70\n",
            nrtk_decimal_year(sc.t0));
    fprintf(fp,"  antex: %s\n",atxfile.c_str());
    fprintf(fp,"stations:\n");
    for (size_t k=0;k<sc.sta.size();k++) {
        auto &s=sc.sta[k];
        fprintf(fp,"  - id: %s\n    ecef: [%.4f, %.4f, %.4f]\n",s.id.c_str(),s.pos[0],s.pos[1],s.pos[2]);
        fprintf(fp,"    antenna: \"%s\"\n    arp_height_m: %.4f\n",s.antenna.c_str(),s.antdel[2]);
        fprintf(fp,"    rinex: %s/%s.rnx\n    rtcm_id: %d\n",sc.outdir.c_str(),s.id.c_str(),(int)(101+k));
    }
    fprintf(fp,"offline:\n  nav: [%s/brdc.nav]\n",sc.outdir.c_str());
    fprintf(fp,"  start: \"%04.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\"\n",ep[0],ep[1],ep[2],ep[3],ep[4],ep[5]);
    fprintf(fp,"  duration_s: %.0f\n  interval_s: %.1f\n",sc.dur,sc.tint);
    fclose(fp);
    fprintf(stderr,"wrote %s (%zu stations)\n",sc.outdir.c_str(),sc.sta.size());
    nrtk_atx_free(&atx);
    return 0;
}
