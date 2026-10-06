/*------------------------------------------------------------------------------
 * ingest.c : RINEX (offline) and RTCM 3 stream (real time) input
 *----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "nrtk/ingest.h"

/* offline RINEX observation source ------------------------------------------*/
extern int nrtk_rnxsrc_open(nrtk_rnxsrc_t *s, const char *file, int rcv)
{
    const char *path=file;
    int stat;

    memset(s,0,sizeof(*s));
    snprintf(s->path,sizeof(s->path),"%s",file);
    s->rcv=rcv;
    if (!init_rnxctr(&s->rnx)) return 0;

    /* compressed file (gzip/Z/Hatanaka): uncompress to a temporary file */
    snprintf(s->tmp,sizeof(s->tmp),"%s",file);
    if ((stat=rtk_uncompress(file,s->tmp))>0) path=s->tmp;
    else {
        s->tmp[0]='\0';
        if (stat<0) {
            nrtk_log(NRTK_LOG_ERROR,"ingest","uncompress failed: %s",file);
            free_rnxctr(&s->rnx);
            return 0;
        }
    }
    if (!(s->fp=fopen(path,"r"))) {
        nrtk_log(NRTK_LOG_ERROR,"ingest","cannot open %s",path);
        free_rnxctr(&s->rnx);
        return 0;
    }
    if (!open_rnxctr(&s->rnx,s->fp)||s->rnx.type!='O') {
        nrtk_log(NRTK_LOG_ERROR,"ingest","not a RINEX observation file: %s",path);
        nrtk_rnxsrc_close(s);
        return 0;
    }
    return 1;
}
extern void nrtk_rnxsrc_close(nrtk_rnxsrc_t *s)
{
    if (s->fp) fclose(s->fp);
    s->fp=NULL;
    free_rnxctr(&s->rnx);
    if (s->tmp[0]) remove(s->tmp);
    s->tmp[0]='\0';
}
static int fill(nrtk_rnxsrc_t *s)
{
    int stat;
    if (s->buffered) return 1;
    if (s->eof||!s->fp) return 0;
    while ((stat=input_rnxctr(&s->rnx,s->fp))!=1) {
        if (stat<0) {s->eof=1; return 0;}  /* end of file */
        if (stat==0&&feof(s->fp)) {s->eof=1; return 0;}
    }
    s->buffered=1;
    return 1;
}
extern gtime_t nrtk_rnxsrc_peek(nrtk_rnxsrc_t *s)
{
    gtime_t t0={0};
    if (!fill(s)) return t0;
    return s->rnx.obs.data[0].time;
}
extern int nrtk_rnxsrc_read(nrtk_rnxsrc_t *s, gtime_t t, obsd_t *obs, int max)
{
    double dt;
    int i,n;

    for (;;) {
        if (!fill(s)) return -1;
        dt=timediff(s->rnx.obs.data[0].time,t);
        if (dt<-DTTOL) {s->buffered=0; continue;}   /* older: skip */
        if (dt>DTTOL) return 0;                      /* newer: keep */
        break;
    }
    n=s->rnx.obs.n<max?s->rnx.obs.n:max;
    for (i=0;i<n;i++) {
        obs[i]=s->rnx.obs.data[i];
        obs[i].rcv=(uint8_t)s->rcv;
    }
    s->buffered=0;
    return n;
}
extern const sta_t *nrtk_rnxsrc_sta(const nrtk_rnxsrc_t *s)
{
    return &s->rnx.sta;
}
/* navigation data -------------------------------------------------------------*/
extern int nrtk_nav_read(const char *files, nav_t *nav)
{
    char buff[4096],*p,*q;
    gtime_t t0={0};
    int nfile=0;

    snprintf(buff,sizeof(buff),"%s",files?files:"");
    for (p=buff;*p;p=q) {
        while (*p==' '||*p==',') p++;
        for (q=p;*q&&*q!=' '&&*q!=',';q++) ;
        if (*q) *q++='\0';
        if (!*p) continue;
        if (readrnxt(p,1,t0,t0,0.0,"",NULL,nav,NULL)<=0) {
            nrtk_log(NRTK_LOG_WARN,"ingest","navigation file read error: %s",p);
            continue;
        }
        nfile++;
    }
    if (nav->n>0) uniqnav(nav);
    nrtk_log(NRTK_LOG_INFO,"ingest","navigation: %d files, %d ephemerides",nfile,nav->n);
    return nav->n>0;
}
extern void nrtk_nav_free(nav_t *nav)
{
    free(nav->eph);  nav->eph=NULL;  nav->n=nav->nmax=0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    free(nav->alm);  nav->alm=NULL;  nav->na=nav->namax=0;
    free(nav->tec);  nav->tec=NULL;  nav->nt=nav->ntmax=0;
}
extern int nrtk_nav_copy(nav_t *dst, const nav_t *src)
{
    eph_t *eph=dst->eph;
    geph_t *geph=dst->geph;
    int nmax=dst->nmax,ngmax=dst->ngmax;

    /* keep destination buffers, copy scalars */
    memcpy(dst,src,sizeof(nav_t));
    dst->eph=eph; dst->nmax=nmax; dst->geph=geph; dst->ngmax=ngmax;
    dst->seph=NULL; dst->ns=dst->nsmax=0;
    dst->peph=NULL; dst->ne=dst->nemax=0;
    dst->pclk=NULL; dst->nc=dst->ncmax=0;
    dst->alm=NULL;  dst->na=dst->namax=0;
    dst->tec=NULL;  dst->nt=dst->ntmax=0;
    if (src->n>dst->nmax) {
        if (!(dst->eph=(eph_t *)realloc(dst->eph,sizeof(eph_t)*src->n))) return 0;
        dst->nmax=src->n;
    }
    if (src->ng>dst->ngmax) {
        if (!(dst->geph=(geph_t *)realloc(dst->geph,sizeof(geph_t)*src->ng))) return 0;
        dst->ngmax=src->ng;
    }
    if (src->n>0) memcpy(dst->eph,src->eph,sizeof(eph_t)*src->n);
    if (src->ng>0) memcpy(dst->geph,src->geph,sizeof(geph_t)*src->ng);
    dst->n=src->n; dst->ng=src->ng;
    return 1;
}
extern int nrtk_nav_add_eph(nav_t *nav, const eph_t *eph)
{
    int i;
    eph_t *p;

    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=eph->sat) continue;
        if (nav->eph[i].iode==eph->iode&&timediff(nav->eph[i].toe,eph->toe)==0.0&&
            nav->eph[i].code==eph->code) {
            nav->eph[i]=*eph;       /* same data set (update health etc.) */
            return 0;
        }
    }
    if (nav->n>=nav->nmax) {
        int nmax=nav->nmax<=0?1024:nav->nmax*2;
        if (!(p=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nmax))) return -1;
        nav->eph=p; nav->nmax=nmax;
    }
    nav->eph[nav->n++]=*eph;
    /* bound memory: drop ephemerides older than 1 day */
    if (nav->n>4000) {
        int j;
        gtime_t tmax=eph->toe;
        for (i=j=0;i<nav->n;i++) {
            if (timediff(tmax,nav->eph[i].toe)<86400.0) nav->eph[j++]=nav->eph[i];
        }
        nav->n=j;
    }
    return 1;
}
/* real-time RTCM 3 stream -----------------------------------------------------*/
extern int nrtk_url2str(const char *url, char *path, int size)
{
    const char *p;

    if (!strncmp(url,"ntrip://",8)) {
        snprintf(path,size,"%s",url+8);
        return STR_NTRIPCLI;
    }
    if (!strncmp(url,"tcpcli://",9)) {snprintf(path,size,"%s",url+9); return STR_TCPCLI;}
    if (!strncmp(url,"tcpsvr://",9)) {snprintf(path,size,"%s",url+9); return STR_TCPSVR;}
    if (!strncmp(url,"serial://",9)) {snprintf(path,size,"%s",url+9); return STR_SERIAL;}
    if (!strncmp(url,"file://",7))   {
        p=url+7;
        snprintf(path,size,"%s",p);
        return STR_FILE;
    }
    return 0;
}
extern int nrtk_strsrc_open(nrtk_strsrc_t *s, const char *url)
{
    char path[1024];
    int type;

    memset(s,0,sizeof(*s));
    snprintf(s->url,sizeof(s->url),"%s",url);
    if (!(type=nrtk_url2str(url,path,sizeof(path)))) {
        nrtk_log(NRTK_LOG_ERROR,"ingest","unsupported stream url: %s",url);
        return 0;
    }
    strinit(&s->str);
    if (!init_rtcm(&s->rtcm)) return 0;
    if (!stropen(&s->str,type,STR_MODE_R,path)) {
        nrtk_log(NRTK_LOG_ERROR,"ingest","stream open error: %s",url);
        free_rtcm(&s->rtcm);
        return 0;
    }
    s->opened=1;
    return 1;
}
extern void nrtk_strsrc_close(nrtk_strsrc_t *s)
{
    if (!s->opened) return;
    strclose(&s->str);
    free_rtcm(&s->rtcm);
    s->opened=0;
}
extern int nrtk_strsrc_poll(nrtk_strsrc_t *s,
                            void (*cb_obs)(void *, const obsd_t *, int),
                            void (*cb_eph)(void *, const eph_t *),
                            void (*cb_sta)(void *, const sta_t *),
                            void *ctx)
{
    uint8_t buff[8192];
    int i,n,ret;

    if (!s->opened) return -1;
    if ((n=strread(&s->str,buff,sizeof(buff)))<=0) {
        strstat(&s->str,s->msg);
        return n;
    }
    s->nbyte+=n;
    for (i=0;i<n;i++) {
        ret=input_rtcm3(&s->rtcm,buff[i]);
        if (ret==1&&s->rtcm.obs.n>0) {
            s->nepoch++;
            s->last_rx=utc2gpst(timeget());
            if (cb_obs) cb_obs(ctx,s->rtcm.obs.data,s->rtcm.obs.n);
        }
        else if (ret==2&&cb_eph) {
            int sat=s->rtcm.ephsat,set=s->rtcm.ephset;
            if (sat>0&&satsys(sat,NULL)!=SYS_GLO) {
                cb_eph(ctx,s->rtcm.nav.eph+(sat-1)+MAXSAT*set);
            }
        }
        else if (ret==5&&cb_sta) {
            cb_sta(ctx,&s->rtcm.sta);
        }
    }
    return n;
}
