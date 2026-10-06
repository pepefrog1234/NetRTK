/*------------------------------------------------------------------------------
 * prep.c : station epoch preparation (signal plan, antenna, clock, geometry)
 *----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>
#include "nrtk/prep.h"

extern void nrtk_prep_opt_default(nrtk_prep_opt_t *opt)
{
    memset(opt,0,sizeof(*opt));
    opt->elmin=5.0*D2R;
    opt->sysmask=(1<<NRTK_GPS)|(1<<NRTK_GAL)|(1<<NRTK_BDS)|(1<<NRTK_QZS);
    opt->tide=1;
    opt->pcv=1;
    nrtk_sigplan_default(opt->plan);
}
/* priority of code in blank separated list (0: not listed, larger is better) */
static int code_pri(uint8_t code, const char *list)
{
    const char *obs=code2obs(code),*p;
    int pri=100;
    size_t n;

    if (!obs||!*obs) return 0;
    n=strlen(obs);
    for (p=list;*p;pri--) {
        while (*p==' '||*p==',') p++;
        if (!*p) break;
        if (!strncmp(p,obs,n)&&(p[n]=='\0'||p[n]==' '||p[n]==',')) return pri;
        while (*p&&*p!=' '&&*p!=',') p++;
    }
    return 0;
}
static void copy_slot(obsd_t *dst, int k, const obsd_t *src, int j)
{
    dst->code[k]=src->code[j];
    dst->L[k]=src->L[j];
    dst->P[k]=src->P[j];
    dst->D[k]=src->D[j];
    dst->SNR[k]=src->SNR[j];
    dst->LLI[k]=src->LLI[j];
    dst->Lstd[k]=src->Lstd[j];
    dst->Pstd[k]=src->Pstd[j];
}
/* reorder signals of one satellite according to the signal plan -------------*/
static int plan_obs(const nrtk_sigplan_t *plan, const obsd_t *raw, obsd_t *out,
                    int *nproc)
{
    int used[NRTK_MAXSIG]={0},k,j,best,pri,bpri,n=0;

    memset(out,0,sizeof(*out));
    out->time=raw->time; out->sat=raw->sat; out->rcv=raw->rcv;
    out->freq=raw->freq;

    for (k=0;k<plan->nproc;k++) {
        for (j=0,best=-1,bpri=0;j<NRTK_MAXSIG;j++) {
            if (used[j]||raw->code[j]==CODE_NONE) continue;
            if (raw->L[j]==0.0&&raw->P[j]==0.0) continue;
            if ((pri=code_pri(raw->code[j],plan->proc[k]))>bpri) {bpri=pri; best=j;}
        }
        if (best>=0) {copy_slot(out,k,raw,best); used[best]=1; n++;}
    }
    *nproc=plan->nproc;
    /* other output signals in the remaining slots */
    for (k=plan->nproc,j=0;j<NRTK_MAXSIG&&k<NRTK_MAXSIG;j++) {
        if (used[j]||raw->code[j]==CODE_NONE) continue;
        if (raw->L[j]==0.0&&raw->P[j]==0.0) continue;
        if (!nrtk_sig_in_list(raw->code[j],plan->out)) continue;
        copy_slot(out,k++,raw,j); used[j]=1; n++;
    }
    return n;
}
static int cmp_sat(const void *a, const void *b)
{
    return (int)((const nrtk_satdata_t *)a)->sat-(int)((const nrtk_satdata_t *)b)->sat;
}
static int cmp_dbl(const void *a, const void *b)
{
    double x=*(const double *)a,y=*(const double *)b;
    return x<y?-1:(x>y?1:0);
}
extern int nrtk_prep_epoch(const nrtk_prep_opt_t *opt, const nrtk_station_t *sta,
                           int ista, const nrtk_ant_t *ant, const obsd_t *obs,
                           int n, const nav_t *nav, nrtk_staepoch_t *ep)
{
    double r0[3],pos[3],e[3],clk[NRTK_MAXSATEP],dtr=0.0,lam;
    gtime_t t,trx;
    int i,j,k,sys,idx,ns=0,nc=0,iter;

    memset(ep,0,sizeof(*ep));
    if (n<=0) return 0;
    t=obs[0].time;
    ep->time=t; ep->sta=ista;

    nrtk_sta_pos(sta,t,r0);
    if (opt->tide) nrtk_tide_pos(t,r0,ep->rr);
    else for (i=0;i<3;i++) ep->rr[i]=r0[i];
    ecef2pos(ep->rr,pos);

    /* signal plan */
    for (i=0;i<n&&ns<NRTK_MAXSATEP;i++) {
        nrtk_satdata_t *s=ep->sat+ns;
        sys=satsys(obs[i].sat,NULL);
        if ((idx=nrtk_sysidx(sys))<0||!(opt->sysmask&(1<<idx))) continue;
        memset(s,0,sizeof(*s));
        s->sat=obs[i].sat; s->sysidx=idx;
        if (plan_obs(&opt->plan[idx],obs+i,&s->obs,&s->nproc)<=0) continue;
        for (k=0;k<NRTK_MAXSIG;k++) s->freq[k]=nrtk_code2freq(sys,s->obs.code[k]);
        ns++;
    }
    /* receiver clock from code with the known position (two passes) */
    for (iter=0;iter<2;iter++) {
        trx=timeadd(t,-dtr);
        for (i=nc=0;i<ns;i++) {
            nrtk_satdata_t *s=ep->sat+i;
            s->r=nrtk_satrange(trx,ep->rr,s->sat,nav,s->rs,s->dts,&s->rate,&s->svh);
            if (s->r<=0.0) continue;
            for (j=0;j<NRTK_MAXSIG;j++) {
                if (s->obs.P[j]!=0.0) {
                    clk[nc++]=(s->obs.P[j]-s->r+CLIGHT*s->dts[0])/CLIGHT;
                    break;
                }
            }
        }
        if (nc<=0) break;
        qsort(clk,nc,sizeof(double),cmp_dbl);
        dtr=nc%2?clk[nc/2]:0.5*(clk[nc/2-1]+clk[nc/2]);
        if (fabs(dtr)<1E-7&&iter==0) break; /* clock steered receiver */
    }
    ep->dtr=dtr;

    /* geometry, health, elevation mask and antenna correction */
    for (i=k=0;i<ns;i++) {
        nrtk_satdata_t s=ep->sat[i];
        if (s.r<=0.0) continue;
        if (satexclude(s.sat,0.0,s.svh,NULL)) continue;
        for (j=0;j<3;j++) e[j]=(s.rs[j]-ep->rr[j])/s.r;
        satazel(pos,e,s.azel);
        if (s.azel[1]<opt->elmin) continue;

        for (j=0;j<NRTK_MAXSIG;j++) {
            if (s.obs.code[j]==CODE_NONE||s.freq[j]<=0.0) continue;
            s.dant[j]=nrtk_ant_dr(ant,nrtk_idxsys(s.sysidx),s.freq[j],s.azel,
                                  sta->antdel,opt->pcv);
            lam=CLIGHT/s.freq[j];
            if (s.obs.L[j]!=0.0) s.obs.L[j]-=s.dant[j]/lam;
            if (s.obs.P[j]!=0.0) s.obs.P[j]-=s.dant[j];
        }
        ep->sat[k++]=s;
    }
    ep->ns=k;
    qsort(ep->sat,k,sizeof(nrtk_satdata_t),cmp_sat);
    return k;
}
extern const nrtk_satdata_t *nrtk_ep_sat(const nrtk_staepoch_t *ep, int sat)
{
    int lo=0,hi=ep->ns-1,mid;

    while (lo<=hi) { /* satellites are sorted by number */
        mid=(lo+hi)/2;
        if (ep->sat[mid].sat==sat) return ep->sat+mid;
        if (ep->sat[mid].sat<sat) lo=mid+1; else hi=mid-1;
    }
    for (mid=0;mid<ep->ns;mid++) if (ep->sat[mid].sat==sat) return ep->sat+mid;
    return NULL;
}
extern double nrtk_phase_m(const nrtk_satdata_t *s, int slot)
{
    if (slot<0||slot>=NRTK_MAXSIG||s->obs.L[slot]==0.0||s->freq[slot]<=0.0) return 0.0;
    return s->obs.L[slot]*CLIGHT/s->freq[slot];
}
extern double nrtk_code_m(const nrtk_satdata_t *s, int slot)
{
    if (slot<0||slot>=NRTK_MAXSIG) return 0.0;
    return s->obs.P[slot];
}
