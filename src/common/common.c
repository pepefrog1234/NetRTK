/*------------------------------------------------------------------------------
 * common.c : libnrtk common functions
 *----------------------------------------------------------------------------*/
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <time.h>
#include "nrtk/common.h"

static const int   sys_tbl[NRTK_NSYS]={SYS_GPS,SYS_GAL,SYS_CMP,SYS_QZS};
static const char *sys_names[NRTK_NSYS]={"GPS","GAL","BDS","QZS"};

extern int nrtk_sysidx(int sys)
{
    switch (sys) {
        case SYS_GPS: return NRTK_GPS;
        case SYS_GAL: return NRTK_GAL;
        case SYS_CMP: return NRTK_BDS;
        case SYS_QZS: return NRTK_QZS;
    }
    return -1;
}
extern int nrtk_idxsys(int idx)
{
    return idx>=0&&idx<NRTK_NSYS?sys_tbl[idx]:SYS_NONE;
}
extern const char *nrtk_sysname(int idx)
{
    return idx>=0&&idx<NRTK_NSYS?sys_names[idx]:"???";
}
extern int nrtk_sysmask(const char *list)
{
    char buff[256],*p,*q;
    int i,mask=0;

    snprintf(buff,sizeof(buff),"%s",list?list:"");
    for (p=buff;*p;p=q) {
        while (*p&&(*p==','||*p==' '||*p=='['||*p==']')) p++;
        for (q=p;*q&&*q!=','&&*q!=' '&&*q!=']';q++) ;
        if (*q) *q++='\0';
        if (!*p) continue;
        for (i=0;i<NRTK_NSYS;i++) {
            if (!strcasecmp(p,sys_names[i])||
                (i==NRTK_BDS&&!strcasecmp(p,"CMP"))||
                (i==NRTK_BDS&&!strcasecmp(p,"BEIDOU"))||
                (i==NRTK_GAL&&!strcasecmp(p,"GALILEO"))||
                (i==NRTK_QZS&&!strcasecmp(p,"QZSS"))) mask|=1<<i;
        }
    }
    return mask;
}
/* default signal plan ---------------------------------------------------------
* processing pairs: GPS/QZS L1-L2 (EWL L2-L5), Galileo E1-E5b (EWL E5b-E5a),
* BeiDou B1I-B3I (EWL B3I-B2a). Output signals follow spec ch.7.
* Galileo uses E5b as second frequency: real data (Leica GR50, Hong Kong
* SatRef 2026-270) showed satellite-dependent half-cycle offsets of E5a phase
* between stations (E5a-E5b MW at 0.3-0.5 cycles, E1-E5b integer).          */
extern void nrtk_sigplan_default(nrtk_sigplan_t plan[NRTK_NSYS])
{
    memset(plan,0,sizeof(nrtk_sigplan_t)*NRTK_NSYS);

    plan[NRTK_GPS].sys=SYS_GPS; plan[NRTK_GPS].nproc=3;
    strcpy(plan[NRTK_GPS].proc[0],"1C 1W 1X 1L");
    strcpy(plan[NRTK_GPS].proc[1],"2W 2L 2X 2S 2P 2C 2D");
    strcpy(plan[NRTK_GPS].proc[2],"5Q 5X 5I");
    strcpy(plan[NRTK_GPS].out,"1C 2W 2L 2X 2S 5Q 5X 5I");
    plan[NRTK_GPS].pfreq[0]=FREQL1; plan[NRTK_GPS].pfreq[1]=FREQL2;
    plan[NRTK_GPS].pfreq[2]=FREQL5;

    plan[NRTK_GAL].sys=SYS_GAL; plan[NRTK_GAL].nproc=3;
    strcpy(plan[NRTK_GAL].proc[0],"1C 1X 1B");
    strcpy(plan[NRTK_GAL].proc[1],"7Q 7X 7I");
    strcpy(plan[NRTK_GAL].proc[2],"5Q 5X 5I");
    strcpy(plan[NRTK_GAL].out,"1C 1X 5Q 5X 7Q 7X");
    plan[NRTK_GAL].pfreq[0]=FREQL1; plan[NRTK_GAL].pfreq[1]=FREQE5b;
    plan[NRTK_GAL].pfreq[2]=FREQL5;

    plan[NRTK_BDS].sys=SYS_CMP; plan[NRTK_BDS].nproc=3;
    strcpy(plan[NRTK_BDS].proc[0],"2I 2X 2Q");
    strcpy(plan[NRTK_BDS].proc[1],"6I 6X 6Q");
    strcpy(plan[NRTK_BDS].proc[2],"5P 5X 5D");
    strcpy(plan[NRTK_BDS].out,"2I 2X 6I 6X 1P 1X 1D 5P 5X 5D");
    plan[NRTK_BDS].pfreq[0]=FREQ1_CMP; plan[NRTK_BDS].pfreq[1]=FREQ3_CMP;
    plan[NRTK_BDS].pfreq[2]=FREQL5;

    plan[NRTK_QZS].sys=SYS_QZS; plan[NRTK_QZS].nproc=3;
    strcpy(plan[NRTK_QZS].proc[0],"1C 1X 1L");
    strcpy(plan[NRTK_QZS].proc[1],"2L 2X 2S");
    strcpy(plan[NRTK_QZS].proc[2],"5Q 5X 5I");
    strcpy(plan[NRTK_QZS].out,"1C 2L 2X 5Q 5X");
    plan[NRTK_QZS].pfreq[0]=FREQL1; plan[NRTK_QZS].pfreq[1]=FREQL2;
    plan[NRTK_QZS].pfreq[2]=FREQL5;
}
extern double nrtk_code2freq(int sys, uint8_t code)
{
    if (code==CODE_NONE) return 0.0;
    return code2freq(sys,code,0);
}
/* test whether signal code is in a blank separated list ("1C 2W ...") ------*/
extern int nrtk_sig_in_list(uint8_t code, const char *list)
{
    const char *obs=code2obs(code),*p;
    size_t n;

    if (!obs||!*obs||!list) return 0;
    n=strlen(obs);
    for (p=list;*p;) {
        while (*p==' '||*p==',') p++;
        if (!strncmp(p,obs,n)&&(p[n]=='\0'||p[n]==' '||p[n]==',')) return 1;
        while (*p&&*p!=' '&&*p!=',') p++;
    }
    return 0;
}
/* decimal year ----------------------------------------------------------------*/
extern double nrtk_decimal_year(gtime_t t)
{
    double ep[6],ep0[6]={0,1,1,0,0,0},ep1[6]={0,1,1,0,0,0};
    gtime_t t0,t1;

    time2epoch(t,ep);
    ep0[0]=ep[0]; ep1[0]=ep[0]+1.0;
    t0=epoch2time(ep0); t1=epoch2time(ep1);
    return ep[0]+timediff(t,t0)/timediff(t1,t0);
}
/* station marker position at time t (velocity applied) -----------------------*/
extern void nrtk_sta_pos(const nrtk_station_t *sta, gtime_t t, double *rr)
{
    double dy=0.0;
    int i;

    if (sta->coord_epoch>0.0&&t.time) dy=nrtk_decimal_year(t)-sta->coord_epoch;
    for (i=0;i<3;i++) rr[i]=sta->pos[i]+sta->vel[i]*dy;
}
/* satellite range by light-time iteration -------------------------------------*/
extern double nrtk_satrange(gtime_t trx, const double *rr, int sat,
                            const nav_t *nav, double *rs, double *dts,
                            double *rate, int *svh)
{
    double r=0.0,tau=0.075,e[3],var=0.0,dr[3],rho;
    gtime_t ttx;
    int i,iter,h=0;

    for (iter=0;iter<6;iter++) {
        ttx=timeadd(trx,-tau);
        if (!satpos(ttx,ttx,sat,EPHOPT_BRDC,nav,rs,dts,&var,&h)) return -1.0;
        for (i=0;i<3;i++) dr[i]=rs[i]-rr[i];
        rho=norm(dr,3);
        if (rho<=0.0) return -1.0;
        if (fabs(rho/CLIGHT-tau)<1E-13) {tau=rho/CLIGHT; break;}
        tau=rho/CLIGHT;
    }
    for (i=0;i<3;i++) e[i]=dr[i]/rho;
    /* Sagnac effect (earth rotation during signal propagation) */
    r=rho+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
    if (rate) {
        *rate=dot3(rs+3,e)+OMGE*(rs[3]*rr[1]-rs[4]*rr[0])/CLIGHT;
    }
    if (svh) *svh=h;
    return r;
}
/* solid earth tide ------------------------------------------------------------*/
extern void nrtk_tide_pos(gtime_t t, const double *r0, double *rr)
{
    double dr[3]={0};
    int i;

    tidedisp(gpst2utc(t),r0,1,NULL,NULL,dr);
    for (i=0;i<3;i++) rr[i]=r0[i]+dr[i];
}
extern void nrtk_ne(const double *r0, const double *r, double *dn, double *de)
{
    double pos[3],d[3],enu[3];
    int i;

    ecef2pos(r0,pos);
    for (i=0;i<3;i++) d[i]=r[i]-r0[i];
    ecef2enu(pos,d,enu);
    *de=enu[0]; *dn=enu[1];
}
extern double nrtk_dist(const double *a, const double *b)
{
    double d[3];
    int i;
    for (i=0;i<3;i++) d[i]=a[i]-b[i];
    return norm(d,3);
}
/* logging ---------------------------------------------------------------------*/
static int log_level=NRTK_LOG_INFO,log_json=0;
static pthread_mutex_t log_lock=PTHREAD_MUTEX_INITIALIZER;
static const char *lvl_str[]={"error","warn","info","debug"};

extern void nrtk_log_level(int level) {log_level=level;}
extern void nrtk_log_json(int on) {log_json=on;}
extern int nrtk_log_enabled(int level) {return level<=log_level;}

extern void nrtk_log(int level, const char *event, const char *fmt, ...)
{
    char msg[2048],ts[64],esc[4096],*q=esc;
    const char *p;
    struct timespec now;
    struct tm tm;
    va_list ap;

    if (level>log_level) return;
    va_start(ap,fmt); vsnprintf(msg,sizeof(msg),fmt,ap); va_end(ap);
    clock_gettime(CLOCK_REALTIME,&now);
    gmtime_r(&now.tv_sec,&tm);
    snprintf(ts,sizeof(ts),"%04d-%02d-%02dT%02d:%02d:%02d.%03dZ",tm.tm_year+1900,
             tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec,
             (int)(now.tv_nsec/1000000));
    pthread_mutex_lock(&log_lock);
    if (log_json) {
        for (p=msg;*p&&q<esc+sizeof(esc)-3;p++) {
            if (*p=='"'||*p=='\\') *q++='\\';
            if (*p=='\n') {*q++='\\'; *q++='n'; continue;}
            *q++=*p;
        }
        *q='\0';
        fprintf(stderr,"{\"ts\":\"%s\",\"level\":\"%s\",\"event\":\"%s\",\"msg\":\"%s\"}\n",
                ts,lvl_str[level<0?0:level>3?3:level],event?event:"",esc);
    }
    else {
        fprintf(stderr,"%s %-5s [%s] %s\n",ts,lvl_str[level<0?0:level>3?3:level],
                event?event:"",msg);
    }
    pthread_mutex_unlock(&log_lock);
}
