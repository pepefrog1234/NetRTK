/*------------------------------------------------------------------------------
 * vrsgen.c : virtual reference station observation generator (spec ch.7)
 *----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>
#include "nrtk/vrsgen.h"

typedef struct {
    int     valid;              /* output in the previous epoch */
    double  I,D;                /* single-difference corrections (m) */
    uint32_t blmask;            /* contributing baselines */
    gtime_t t;
} prevsat_t;

struct nrtk_vrs {
    nrtk_vrs_opt_t opt;
    int     master;             /* current master / pass-through station */
    int     mode;               /* 0: none, 1: VRS, 2: pass-through */
    int     refsat[NRTK_NSYS];
    double  offI[NRTK_NSYS],offD[NRTK_NSYS];
    prevsat_t prev[MAXSAT];
};

extern void nrtk_vrs_opt_default(nrtk_vrs_opt_t *o)
{
    memset(o,0,sizeof(*o));
    o->elmin=10.0*D2R;
    o->sysmask=(1<<NRTK_GPS)|(1<<NRTK_GAL)|(1<<NRTK_BDS)|(1<<NRTK_QZS);
    o->pass_unfixed=0;
    o->jump_thres=0.03;
    o->tide=1;
    nrtk_sigplan_default(o->plan);
}
extern nrtk_vrs_t *nrtk_vrs_create(const nrtk_vrs_opt_t *opt)
{
    nrtk_vrs_t *v=(nrtk_vrs_t *)calloc(1,sizeof(nrtk_vrs_t));
    if (!v) return NULL;
    v->opt=*opt;
    v->master=-1;
    return v;
}
extern void nrtk_vrs_free(nrtk_vrs_t *v) {free(v);}

extern void nrtk_vrs_reset(nrtk_vrs_t *v)
{
    v->master=-1; v->mode=0;
    memset(v->refsat,0,sizeof(v->refsat));
    memset(v->offI,0,sizeof(v->offI));
    memset(v->offD,0,sizeof(v->offD));
    memset(v->prev,0,sizeof(v->prev));
}
extern const nrtk_vrs_opt_t *nrtk_vrs_getopt(const nrtk_vrs_t *v) {return &v->opt;}

static const nrtk_vsat_t *corr_sat(const nrtk_vcorr_t *c, int sat)
{
    int i;
    for (i=0;i<c->ns;i++) if (c->sat[i].sat==sat) return c->sat+i;
    return NULL;
}
/* single-difference offsets that keep the corrections continuous when the
 * DD reference changes: off = mean(prevSD(s) - DD(s)) over continued sats */
static void update_offsets(nrtk_vrs_t *v, const nrtk_vcorr_t *c)
{
    int sys,i,n;
    double sI,sD;

    for (sys=0;sys<NRTK_NSYS;sys++) {
        if (!c->refsat[sys]) continue;
        if (v->refsat[sys]==c->refsat[sys]) continue;
        for (i=n=0,sI=sD=0.0;i<c->ns;i++) {
            const nrtk_vsat_t *s=c->sat+i;
            const prevsat_t *p=v->prev+s->sat-1;
            if (s->sysidx!=sys||!s->ok||!p->valid) continue;
            sI+=p->I-s->dd_iono; sD+=p->D-s->dd_nondisp; n++;
        }
        v->offI[sys]=n>0?sI/n:0.0;
        v->offD[sys]=n>0?sD/n:0.0;
        v->refsat[sys]=c->refsat[sys];
    }
}
static int out_sig(const nrtk_vrs_opt_t *o, int sysidx, uint8_t code)
{
    return code!=CODE_NONE&&nrtk_sig_in_list(code,o->plan[sysidx].out);
}
extern int nrtk_vrs_generate(nrtk_vrs_t *v, const nrtk_staepoch_t *mep,
                             const nrtk_vcorr_t *corr, const nav_t *nav,
                             obsd_t *out, int maxn, nrtk_vrs_info_t *info)
{
    double rv[3],posv[3],rs[6],dts[2],e[3],azel[2],rhov,ratev,lam,mu,dI,dD,dr;
    gtime_t trx;
    int i,j,k,n=0,nreset=0,reset,svh,newmaster;
    prevsat_t cur[MAXSAT];

    memset(cur,0,sizeof(cur));
    if (info) memset(info,0,sizeof(*info));
    newmaster=v->mode!=1||v->master!=mep->sta;
    if (newmaster) {
        nrtk_vrs_reset(v);
        v->master=mep->sta; v->mode=1;
    }
    update_offsets(v,corr);

    if (v->opt.tide) nrtk_tide_pos(mep->time,v->opt.pos,rv);
    else memcpy(rv,v->opt.pos,sizeof(rv));
    ecef2pos(rv,posv);
    trx=timeadd(mep->time,-mep->dtr);   /* true reception time of the master */

    for (i=0;i<mep->ns&&n<maxn;i++) {
        const nrtk_satdata_t *s=mep->sat+i;
        const nrtk_vsat_t *c;
        obsd_t *o;
        int corrected;

        if (!(v->opt.sysmask&(1<<s->sysidx))) continue;
        c=corr_sat(corr,s->sat);
        corrected=c&&c->ok;
        if (!corrected&&!v->opt.pass_unfixed) continue;

        if ((rhov=nrtk_satrange(trx,rv,s->sat,nav,rs,dts,&ratev,&svh))<=0.0) continue;
        for (k=0;k<3;k++) e[k]=(rs[k]-rv[k])/rhov;
        satazel(posv,e,azel);
        if (azel[1]<v->opt.elmin) continue;

        dr=rhov-s->r;
        dI=dD=0.0;
        if (corrected) {
            dI=c->dd_iono+v->offI[s->sysidx];
            dD=c->dd_nondisp+v->offD[s->sysidx];
        }
        /* lock time reset conditions */
        reset=0;
        if (!corrected) reset=1;
        else {
            const prevsat_t *p=v->prev+s->sat-1;
            if (!p->valid) reset=1;                       /* new satellite / new VRS */
            else if (c->newfix) reset=1;                  /* network (re)fix */
            else if (c->blmask!=p->blmask) reset=1;       /* interpolation basis changed */
            else {
                /* jump of the correction actually applied to each processing
                   frequency (dD - mu_j dI). Testing dD/dI themselves would react
                   to the noise of the ionosphere-free-like split. */
                for (j=0;j<s->nproc&&j<NRTK_NPROC&&!reset;j++) {
                    if (s->freq[j]<=0.0) continue;
                    mu=(FREQL1/s->freq[j])*(FREQL1/s->freq[j]);
                    if (fabs((dD-mu*dI)-(p->D-mu*p->I))>v->opt.jump_thres) reset=1;
                }
            }
        }
        o=out+n;
        memset(o,0,sizeof(*o));
        o->time=mep->time; o->sat=s->sat; o->rcv=1;
        for (j=k=0;j<NRTK_MAXSIG;j++) {
            uint8_t code=s->obs.code[j];
            if (!out_sig(&v->opt,s->sysidx,code)||s->freq[j]<=0.0) continue;
            lam=CLIGHT/s->freq[j];
            mu=(FREQL1/s->freq[j])*(FREQL1/s->freq[j]);
            o->code[j]=code;
            if (s->obs.L[j]!=0.0) o->L[j]=s->obs.L[j]+(dr+dD-mu*dI)/lam;
            if (s->obs.P[j]!=0.0) o->P[j]=s->obs.P[j]+dr+dD+mu*dI;
            if (s->obs.D[j]!=0.0f) o->D[j]=(float)(s->obs.D[j]-(ratev-s->rate)/lam);
            o->SNR[j]=s->obs.SNR[j];
            o->LLI[j]=(uint8_t)(s->obs.LLI[j]|(reset?LLI_SLIP:0));
            o->Lstd[j]=s->obs.Lstd[j]; o->Pstd[j]=s->obs.Pstd[j];
            k++;
        }
        if (k==0) continue;
        if (reset) nreset++;
        if (corrected) {
            cur[s->sat-1].valid=1;
            cur[s->sat-1].I=dI; cur[s->sat-1].D=dD;
            cur[s->sat-1].blmask=c->blmask;
            cur[s->sat-1].t=mep->time;
        }
        n++;
    }
    memcpy(v->prev,cur,sizeof(cur));
    if (info) {
        info->nsat=n; info->nlockreset=nreset; info->master=mep->sta; info->dtr=mep->dtr;
    }
    return n;
}
extern int nrtk_vrs_passthrough(nrtk_vrs_t *v, const nrtk_staepoch_t *sep,
                                obsd_t *out, int maxn, nrtk_vrs_info_t *info)
{
    int i,j,n=0,reset;

    if (info) memset(info,0,sizeof(*info));
    reset=v->mode!=2||v->master!=sep->sta;
    if (reset) {nrtk_vrs_reset(v); v->master=sep->sta; v->mode=2;}

    for (i=0;i<sep->ns&&n<maxn;i++) {
        const nrtk_satdata_t *s=sep->sat+i;
        obsd_t *o=out+n;
        int k=0;
        if (!(v->opt.sysmask&(1<<s->sysidx))||s->azel[1]<v->opt.elmin) continue;
        memset(o,0,sizeof(*o));
        o->time=sep->time; o->sat=s->sat; o->rcv=1;
        for (j=0;j<NRTK_MAXSIG;j++) {
            if (!out_sig(&v->opt,s->sysidx,s->obs.code[j])) continue;
            o->code[j]=s->obs.code[j];
            o->L[j]=s->obs.L[j]; o->P[j]=s->obs.P[j]; o->D[j]=s->obs.D[j];
            o->SNR[j]=s->obs.SNR[j];
            o->LLI[j]=(uint8_t)(s->obs.LLI[j]|(reset||!v->prev[s->sat-1].valid?LLI_SLIP:0));
            k++;
        }
        if (!k) continue;
        n++;
    }
    for (i=0;i<MAXSAT;i++) v->prev[i].valid=0;
    for (i=0;i<n;i++) v->prev[out[i].sat-1].valid=1;
    if (info) {info->nsat=n; info->master=sep->sta; info->nlockreset=reset?n:0; info->dtr=sep->dtr;}
    return n;
}
