/*------------------------------------------------------------------------------
 * netpos.c : network solution between reference stations (spec ch.5)
 *
 * Shared layers (both engines):
 *   - slip detection by MW / geometry-free / LLI / signal change
 *   - extra-wide-lane (f2-f3) and wide-lane (f1-f2) fixing from time-averaged
 *     Melbourne-Wubbena combinations, referred to all fixed satellites
 *   - L1 ambiguity resolution with LAMBDA (ratio test, float std gate,
 *     partial fixing dropping the lowest elevations first)
 *   - validation of new fixes and monitoring of held ones against an
 *     independent troposphere model (Saastamoinen+GMF plus one robust relative
 *     zenith wet delay); a wide-lane error shows as lam_NL*f2/(f1-f2)
 *   - fix-and-hold, ambiguity-fixed DD phase residuals split into the
 *     dispersive (ionosphere, L1) and non-dispersive part
 * Float engines (solver.engine):
 *   rtklib : RTKLIB-EX rtkpos() in fixed-position mode (IONOOPT_EST, ZWD
 *            estimated, ARMODE_OFF), unmodified; slips injected via LLI,
 *            master ZWD pinned to the model (only the relative delay is
 *            observable), integer search and hold done here
 *   native : ionosphere-free phase with the fixed wide-lane removed; states
 *            are one relative ZWD and one float L1 ambiguity per satellite.
 *            Code and ionosphere do not enter the L1 float, so code multipath
 *            and ionospheric activity cannot bias it.
 * Integer ambiguities are kept as single-difference integers with an
 * arbitrary per-system datum: the DD reference satellite can change freely
 * (DD = difference of SD integers) without touching the filters.
 *----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>
#include "nrtk/netpos.h"
#include "nrtk/tropo.h"

/* RTKLIB-EX v2.5.1 rtkpos.c state vector layout (see rtkpos.c NP..IB) */
#define NP_(o)     ((o)->dynamics==0?3:9)
#define NI_(o)     ((o)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT_(o)     ((o)->tropopt<TROPOPT_EST?0:((o)->tropopt<TROPOPT_ESTG?2:6))
#define NL_(o)     ((o)->glomodear!=GLO_ARMODE_AUTOCAL?0:NFREQGLO)
#define NF_(o)     ((o)->ionoopt==IONOOPT_IFLC?1:(o)->nf)
#define NR_(o)     (NP_(o)+NI_(o)+NT_(o)+NL_(o))
#define NX_(o)     (NR_(o)+MAXSAT*NF_(o))
#define IT_(r,o)   (NP_(o)+NI_(o)+NT_(o)/2*(r))  /* r: 0=rover(aux),1=base(master) */
#define IB_(s,f,o) (NR_(o)+MAXSAT*(f)+(s)-1)

#define NAT_NX     (MAXSAT+1)   /* native: 0 = relative ZWD, s = L1 ambiguity of sat s */

#define SQR(x)     ((x)*(x))
#define MIN(x,y)   ((x)<=(y)?(x):(y))
#define MAX(x,y)   ((x)>=(y)?(x):(y))
#define ROUND(x)   floor((x)+0.5)
#define MAXGAP     10.0         /* max data gap keeping ambiguities (s) */
#define VAR_COND   1E-6         /* conditioning constraint variance (cyc^2) */
#define NL_MAXSTD  0.15         /* max float std of a DD L1 ambiguity to fix (cyc) */
#define PAR_AMBIG  0.25         /* partial fixing: drop a DD this far from its integer first (cyc) */
#define REL_ZWD    0.05         /* a-priori std of the aux ZWD w.r.t. the model (m) */
#define DZWD_PRIOR 0.05         /* prior std of the relative ZWD in validation (m) */
#define DZWD_MAX   0.12         /* implausible relative ZWD (m) */
#define JUMP_RES   0.07         /* max epoch-to-epoch change of a DD residual (m) */
#define VAR_NEWAMB 1E4          /* initial variance of a native L1 ambiguity (cyc^2) */
#define HELD_TAU   120.0        /* averaging time of the held-fix model residual (s) */
#define HELD_MIN   60.0         /* min averaging time before the mean is tested (s) */
#define HELD_EL    (25.0*D2R)   /* min elevation of the mean residual test */
#define HELD_THR   0.065        /* max |mean residual| of a held fix (m), ~0.6 lam_NL */

typedef struct {
    int     lock;               /* consecutive epochs with valid data */
    gtime_t tlast;              /* last epoch with valid data */
    uint8_t code[2][NRTK_NPROC];/* signal codes {aux,master} per proc slot */
    double  mw[2],mw2[2];       /* sum/sum^2 of SD MW: [0] f1-f2, [1] f2-f3 (cyc) */
    int     nmw[2];
    gtime_t mwt0[2];            /* start of the MW average */
    int     nrej;               /* consecutive rejected L1 fixes */
    int     nbad;               /* consecutive epochs failing the model check */
    double  gf;                 /* last SD geometry-free phase (m) */
    int     wl,ewl,nl;          /* fixed flags (nl==2: tentative) */
    double  Nwl,Newl,N1;        /* SD integers with per-system datum */
    double  Nwl_nat;            /* wide-lane used by the native float state */
    int     nfixep;
    double  res_prev[3];        /* last DD residual per proc slot (continuity check) */
    int     res_ref;            /* reference sat of res_prev */
    int     res_has3;           /* res_prev[2] valid */
    double  eavg,eavg_w;        /* mean model residual of a held fix (m), averaging time (s) */
    gtime_t eavg_t;
    int     eavg_anc;           /* anchor satellite of eavg */
    int     newfix;
} satst_t;

typedef struct {                /* native ionosphere-free filter */
    double *x,*P;               /* NAT_NX states */
    gtime_t t;
    int     init;
} natf_t;

struct nrtk_baseline {
    nrtk_netpos_opt_t opt;
    int     ista_m,ista_a;
    double  rm[3],ra[3];        /* marker positions (no tide) */
    double  length;
    int     nf;
    int     engine;
    rtk_t   rtk;                /* rtklib engine */
    int     trop_init;
    natf_t  nat;                /* native engine */
    satst_t ss[MAXSAT];
    double  ratio;
    double  dzwd;               /* relative ZWD aux-master w.r.t. the model (m) */
    int     dzwd_ok;
    /* statistics */
    long    nep,nfixep,nslip,nrefix;
    gtime_t tfirst;
    double  ttff;
    int     fixed_prev,everfixed;
    uint8_t hist[60];
    int     hn,hi;
};

extern void nrtk_netpos_opt_default(nrtk_netpos_opt_t *o)
{
    memset(o,0,sizeof(*o));
    o->engine=NRTK_ENGINE_NATIVE;
    o->sysmask=(1<<NRTK_GPS)|(1<<NRTK_GAL)|(1<<NRTK_BDS)|(1<<NRTK_QZS);
    o->nfreq=3;
    o->elmin=10.0*D2R;
    o->elmin_ar=10.0*D2R;
    o->ar_ratio=3.0;
    o->min_fix=4;
    o->min_lock=5;
    o->max_drop=6;
    o->iono_prior_ppm=0.5;
    o->iono_prior_min=1.00;
    o->iono_prn=1E-3;
    o->trop_prn=1E-4;
    o->trop_std=0.15;
    o->trop_grad=0;
    o->err_phase=0.003;
    o->code_ratio=100.0;
    o->wl_min_epochs=5;
    o->wl_min_time=30.0;
    o->wl_max_frac=0.20;
    o->wl_max_std=0.10;
    o->mw_slip=1.0;
    o->gf_slip=0.05;
    o->hold_var=0.001;
    o->max_res=0.06;
    o->tide=1;
    o->closure=1;
}
extern void nrtk_split_iono(double fa, double fb, double da, double db,
                            double *iono_l1, double *nondisp)
{
    double ga=SQR(FREQL1/fa),gb=SQR(FREQL1/fb);
    /* da = D - ga*I, db = D - gb*I */
    double I=(da-db)/(gb-ga);
    *iono_l1=I;
    *nondisp=da+ga*I;
}
/* create ----------------------------------------------------------------------*/
static void setup_prcopt(nrtk_baseline_t *bl, prcopt_t *p)
{
    const nrtk_netpos_opt_t *o=&bl->opt;
    double stdi;
    int i;

    *p=prcopt_default;
    p->mode=PMODE_FIXED;
    p->soltype=0;
    p->nf=bl->nf;
    p->navsys=0;
    for (i=0;i<NRTK_NSYS;i++) if (o->sysmask&(1<<i)) p->navsys|=nrtk_idxsys(i);
    p->elmin=o->elmin;
    p->sateph=EPHOPT_BRDC;
    p->modear=ARMODE_OFF;
    p->glomodear=GLO_ARMODE_OFF;
    p->bdsmodear=1;
    p->arfilter=0;
    p->maxout=5;
    p->minlock=0;
    p->ionoopt=IONOOPT_EST;
    p->tropopt=o->trop_grad?TROPOPT_ESTG:TROPOPT_EST;
    p->dynamics=0;
    p->tidecorr=o->tide?1:0;
    p->niter=1;
    p->refpos=POSOPT_POS_XYZ;
    p->rovpos=POSOPT_POS_XYZ;
    for (i=0;i<NFREQ;i++) p->eratio[i]=o->code_ratio;
    p->err[1]=o->err_phase;
    p->err[2]=o->err_phase;
    p->err[3]=0.0;
    p->err[6]=0.0;
    p->err[7]=0.0;
    /* SD ionosphere prior: RTKLIB uses std[1]*bl/1E4 */
    stdi=MAX(o->iono_prior_ppm*1E-6*bl->length,o->iono_prior_min);
    p->std[0]=30.0;
    p->std[1]=bl->length>1.0?stdi*1E4/bl->length:stdi;
    p->std[2]=o->trop_std;
    p->prn[0]=1E-4;
    p->prn[1]=o->iono_prn;
    p->prn[2]=o->trop_prn;
    p->thresslip=o->gf_slip;
    p->thresdop=0.0;
    p->maxinno[0]=5.0;
    p->maxinno[1]=30.0;
    p->maxtdiff=30.0;
    p->outsingle=0;
    p->posopt[1]=0;             /* observations are antenna corrected */
    memset(p->pcvr,0,sizeof(p->pcvr));
    memset(p->antdel,0,sizeof(p->antdel));
    for (i=0;i<3;i++) {p->ru[i]=bl->ra[i]; p->rb[i]=bl->rm[i];}
}
extern nrtk_baseline_t *nrtk_bl_create(const nrtk_netpos_opt_t *opt,
                                       const nrtk_station_t *m, int ista_m,
                                       const nrtk_station_t *a, int ista_a)
{
    nrtk_baseline_t *bl;
    prcopt_t popt;

    if (!(bl=(nrtk_baseline_t *)calloc(1,sizeof(nrtk_baseline_t)))) return NULL;
    bl->opt=*opt;
    bl->ista_m=ista_m; bl->ista_a=ista_a;
    memcpy(bl->rm,m->pos,sizeof(bl->rm));
    memcpy(bl->ra,a->pos,sizeof(bl->ra));
    bl->length=nrtk_dist(bl->rm,bl->ra);
    bl->nf=opt->nfreq>=3?3:2;
    if (bl->nf>NFREQ) bl->nf=NFREQ;
    bl->engine=opt->engine;
    bl->ttff=-1.0;
    if (bl->engine==NRTK_ENGINE_RTKLIB) {
        setup_prcopt(bl,&popt);
        rtkinit(&bl->rtk,&popt);
        if (bl->rtk.nx!=NX_(&bl->rtk.opt)) {
            nrtk_log(NRTK_LOG_ERROR,"netpos","RTKLIB state layout mismatch (nx=%d expected %d)",
                     bl->rtk.nx,NX_(&bl->rtk.opt));
            rtkfree(&bl->rtk); free(bl);
            return NULL;
        }
    }
    else {
        bl->nat.x=zeros(NAT_NX,1);
        bl->nat.P=zeros(NAT_NX,NAT_NX);
        if (!bl->nat.x||!bl->nat.P) {free(bl->nat.x); free(bl->nat.P); free(bl); return NULL;}
    }
    nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d len=%.1f km nf=%d engine=%s",ista_m,ista_a,
             bl->length/1E3,bl->nf,bl->engine==NRTK_ENGINE_RTKLIB?"rtklib":"native");
    if (opt->trace>0) {traceopen("netpos.trace"); tracelevel(opt->trace);}
    return bl;
}
extern void nrtk_bl_free(nrtk_baseline_t *bl)
{
    if (!bl) return;
    if (bl->engine==NRTK_ENGINE_RTKLIB) rtkfree(&bl->rtk);
    free(bl->nat.x); free(bl->nat.P);
    free(bl);
}
extern void nrtk_bl_reset(nrtk_baseline_t *bl)
{
    if (bl->engine==NRTK_ENGINE_RTKLIB) {
        prcopt_t popt=bl->rtk.opt;
        rtkfree(&bl->rtk);
        rtkinit(&bl->rtk,&popt);
    }
    else {
        memset(bl->nat.x,0,sizeof(double)*NAT_NX);
        memset(bl->nat.P,0,sizeof(double)*NAT_NX*NAT_NX);
        bl->nat.init=0; bl->nat.t.time=0;
    }
    memset(bl->ss,0,sizeof(bl->ss));
    bl->trop_init=0;
    bl->dzwd_ok=0;
}
/* observation combinations ------------------------------------------------------*/
static int has_slots(const nrtk_satdata_t *s, int a, int b)
{
    return s->nproc>b&&s->obs.L[a]!=0.0&&s->obs.L[b]!=0.0&&s->obs.P[a]!=0.0&&
           s->obs.P[b]!=0.0&&s->freq[a]>0.0&&s->freq[b]>0.0;
}
/* Melbourne-Wubbena (cycles of the wide-lane) */
static double mw_cyc(const nrtk_satdata_t *s, int a, int b)
{
    double fa=s->freq[a],fb=s->freq[b],La=nrtk_phase_m(s,a),Lb=nrtk_phase_m(s,b);
    double mw=(fa*La-fb*Lb)/(fa-fb)-(fa*s->obs.P[a]+fb*s->obs.P[b])/(fa+fb);
    return mw/(CLIGHT/(fa-fb));
}
static double gf_m(const nrtk_satdata_t *s)
{
    return nrtk_phase_m(s,0)-nrtk_phase_m(s,1);
}
static void reset_mw(satst_t *st, int k)
{
    st->mw[k]=st->mw2[k]=0.0;
    st->nmw[k]=0;
    st->mwt0[k].time=0; st->mwt0[k].sec=0.0;
}
static void reset_sat(satst_t *st)
{
    st->lock=0;
    reset_mw(st,0); reset_mw(st,1);
    st->nrej=st->nbad=0;
    st->wl=st->ewl=st->nl=0;
    st->nfixep=0;
    st->gf=0.0;
    st->res_ref=0;
    st->eavg_w=0.0;
}
static double mw_mean(const satst_t *st, int k)
{
    return st->nmw[k]>0?st->mw[k]/st->nmw[k]:0.0;
}
static double mw_var_mean(const satst_t *st, int k) /* variance of the mean */
{
    double m,v;
    if (st->nmw[k]<2) return 1.0;
    m=st->mw[k]/st->nmw[k];
    v=(st->mw2[k]-st->nmw[k]*m*m)/(st->nmw[k]-1);
    if (v<0.01) v=0.01;         /* floor: multipath is time correlated */
    return v/st->nmw[k]*5.0;    /* inflate for correlated noise */
}
/* common satellites of the two stations -----------------------------------------*/
typedef struct {
    int sat,sysidx;
    const nrtk_satdata_t *m,*a;
    double el;                  /* elevation at master */
    int ok2,ok3;                /* f1/f2 and f2/f3 available at both */
} comsat_t;

static int common_sats(const nrtk_baseline_t *bl, const nrtk_staepoch_t *em,
                       const nrtk_staepoch_t *ea, comsat_t *cs)
{
    int i,n=0;
    const nrtk_satdata_t *sa;

    for (i=0;i<em->ns;i++) {
        const nrtk_satdata_t *sm=em->sat+i;
        if (!(bl->opt.sysmask&(1<<sm->sysidx))) continue;
        if (!(sa=nrtk_ep_sat(ea,sm->sat))) continue;
        if (sm->azel[1]<bl->opt.elmin||sa->azel[1]<bl->opt.elmin) continue;
        cs[n].sat=sm->sat; cs[n].sysidx=sm->sysidx;
        cs[n].m=sm; cs[n].a=sa; cs[n].el=sm->azel[1];
        cs[n].ok2=has_slots(sm,0,1)&&has_slots(sa,0,1)&&
                  sm->freq[0]==sa->freq[0]&&sm->freq[1]==sa->freq[1];
        cs[n].ok3=bl->nf>=3&&has_slots(sm,1,2)&&has_slots(sa,1,2)&&
                  sm->freq[2]==sa->freq[2];
        n++;
    }
    return n;
}
/* slip detection and MW accumulation ----------------------------------------------*/
static void detect_slips(nrtk_baseline_t *bl, gtime_t t, comsat_t *cs, int n,
                         int *slip)
{
    int i,j,k,lli;
    double dt,gf,mw,lim;

    for (i=0;i<n;i++) {
        satst_t *st=bl->ss+cs[i].sat-1;
        slip[i]=0;
        if (!cs[i].ok2) continue;
        dt=st->tlast.time?timediff(t,st->tlast):0.0;

        lli=0;
        for (j=0;j<bl->nf&&j<NRTK_NPROC;j++) {
            lli|=cs[i].m->obs.LLI[j]|cs[i].a->obs.LLI[j];
            if (st->lock>0&&(st->code[0][j]!=cs[i].a->obs.code[j]||
                             st->code[1][j]!=cs[i].m->obs.code[j])) lli|=LLI_SLIP;
        }
        if (lli&LLI_SLIP) slip[i]=1;
        if (st->lock>0&&dt>MAXGAP) slip[i]=1;

        gf=gf_m(cs[i].a)-gf_m(cs[i].m);
        if (st->lock>0&&dt>0.0&&bl->opt.gf_slip>0.0) {
            lim=bl->opt.gf_slip*(1.0+dt/30.0);
            if (fabs(gf-st->gf)>lim) slip[i]=1;
        }
        for (k=0;k<2;k++) {
            if (k==1&&!cs[i].ok3) continue;
            mw=k==0?mw_cyc(cs[i].a,0,1)-mw_cyc(cs[i].m,0,1):
                    mw_cyc(cs[i].a,1,2)-mw_cyc(cs[i].m,1,2);
            if (st->nmw[k]>=5&&fabs(mw-mw_mean(st,k))>bl->opt.mw_slip*(k?0.5:1.0)+
                4.0*sqrt(mw_var_mean(st,k)*st->nmw[k]/5.0)) {
                slip[i]=1;
            }
        }
        if (slip[i]&&st->lock>0) {
            if (st->nl) {
                bl->nrefix++;
                nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: slip",bl->ista_m,
                         bl->ista_a,cs[i].sat);
            }
            reset_sat(st);
            bl->nslip++;
        }
        /* accumulate */
        st->gf=gf;
        for (k=0;k<2;k++) {
            if (k==1&&!cs[i].ok3) continue;
            mw=k==0?mw_cyc(cs[i].a,0,1)-mw_cyc(cs[i].m,0,1):
                    mw_cyc(cs[i].a,1,2)-mw_cyc(cs[i].m,1,2);
            if (st->nmw[k]==0) st->mwt0[k]=t;
            st->mw[k]+=mw; st->mw2[k]+=mw*mw; st->nmw[k]++;
        }
        for (j=0;j<bl->nf&&j<NRTK_NPROC;j++) {
            st->code[0][j]=cs[i].a->obs.code[j];
            st->code[1][j]=cs[i].m->obs.code[j];
        }
        st->lock++;
        st->tlast=t;
    }
}
/* engine independent access to the float L1 ambiguities ---------------------*/
static int rtk_bias_ok(const rtk_t *rtk, int sat, int f)
{
    int j=IB_(sat,f,&rtk->opt);
    return f<rtk->opt.nf&&rtk->x[j]!=0.0&&rtk->P[j+j*rtk->nx]>0.0;
}
static int amb_ok(const nrtk_baseline_t *bl, int sat)
{
    if (bl->engine==NRTK_ENGINE_RTKLIB) {
        return rtk_bias_ok(&bl->rtk,sat,0)&&rtk_bias_ok(&bl->rtk,sat,1);
    }
    return bl->nat.x[sat]!=0.0&&bl->nat.P[sat*(NAT_NX+1)]>0.0;
}
static int amb_idx(const nrtk_baseline_t *bl, int sat)
{
    return bl->engine==NRTK_ENGINE_RTKLIB?IB_(sat,0,&bl->rtk.opt):sat;
}
static int amb_sat(const nrtk_baseline_t *bl, int idx)
{
    return bl->engine==NRTK_ENGINE_RTKLIB?idx-IB_(1,0,&bl->rtk.opt)+1:idx;
}
static int eng_nx(const nrtk_baseline_t *bl)
{
    return bl->engine==NRTK_ENGINE_RTKLIB?bl->rtk.nx:NAT_NX;
}
static double *eng_x(nrtk_baseline_t *bl)
{
    return bl->engine==NRTK_ENGINE_RTKLIB?bl->rtk.x:bl->nat.x;
}
static double *eng_P(nrtk_baseline_t *bl)
{
    return bl->engine==NRTK_ENGINE_RTKLIB?bl->rtk.P:bl->nat.P;
}
/* drop the float ambiguity of a satellite whose integer was wrong: fix-and-hold
 * has pinned it to that integer, so a refix would just repeat it */
static void nat_clear(natf_t *f, int i);
static void reset_amb(nrtk_baseline_t *bl, int sat)
{
    int f,j,k,nx;
    bl->ss[sat-1].eavg_w=0.0;
    if (bl->engine!=NRTK_ENGINE_RTKLIB) {nat_clear(&bl->nat,sat); return;}
    nx=bl->rtk.nx;
    for (f=0;f<bl->rtk.opt.nf;f++) {      /* rtkpos reinitializes zero biases */
        j=IB_(sat,f,&bl->rtk.opt);
        bl->rtk.x[j]=0.0;
        for (k=0;k<nx;k++) bl->rtk.P[j+k*nx]=bl->rtk.P[k+j*nx]=0.0;
    }
}
/* best (highest elevation) satellite of a system satisfying a flag ----------
 * type 2 (held L1 reference) also needs both phases this epoch: every DD
 * residual of the system is formed against it */
static int anchor(const nrtk_baseline_t *bl, const comsat_t *cs, int n, int sys,
                  int type)
{
    int i,best=-1,rtk=bl->engine==NRTK_ENGINE_RTKLIB;
    for (i=0;i<n;i++) {
        const satst_t *st=bl->ss+cs[i].sat-1;
        if (cs[i].sysidx!=sys) continue;
        if (type==0&&!(st->wl&&cs[i].ok2&&(!rtk||(rtk_bias_ok(&bl->rtk,cs[i].sat,0)&&
                                                 rtk_bias_ok(&bl->rtk,cs[i].sat,1))))) continue;
        if (type==1&&!(st->ewl&&cs[i].ok3&&(!rtk||(rtk_bias_ok(&bl->rtk,cs[i].sat,1)&&
                                                  rtk_bias_ok(&bl->rtk,cs[i].sat,2))))) continue;
        if (type==2&&!(st->nl==1&&cs[i].ok2&&amb_ok(bl,cs[i].sat))) continue;
        if (best<0||cs[i].el>cs[best].el) best=i;
    }
    return best;
}
static double N2of(const satst_t *st) {return st->N1-st->Nwl;}
static double N3of(const satst_t *st) {return st->N1-st->Nwl-st->Newl;}

/* integer constraints ------------------------------------------------------------*/
typedef struct { int i[4]; double h[4]; int n; double v,var; } cons_t;

/* L1 (rtklib: L1/L2/L3) integers of held satellites; with_wl (rtklib only):
 * also the fixed (extra-)wide-lanes, used to condition the ambiguity search */
static int collect_cons(nrtk_baseline_t *bl, const comsat_t *cs, int n,
                        double var_wl, double var_nl, int with_wl, cons_t *c, int maxc)
{
    const double *x=eng_x(bl);
    int sys,i,k,nc=0,w,e,a,nk;
    int rtk=bl->engine==NRTK_ENGINE_RTKLIB;
    const prcopt_t *o=&bl->rtk.opt;

    for (sys=0;sys<NRTK_NSYS;sys++) {
        a=anchor(bl,cs,n,sys,2);
        w=rtk&&with_wl?anchor(bl,cs,n,sys,0):-1;
        e=rtk&&with_wl?anchor(bl,cs,n,sys,1):-1;
        for (i=0;i<n&&nc+5<=maxc;i++) {   /* up to 5 constraints per satellite */
            const satst_t *st=bl->ss+cs[i].sat-1;
            int s=cs[i].sat;
            if (cs[i].sysidx!=sys) continue;
            if (a>=0&&i!=a&&st->nl==1&&amb_ok(bl,s)) {
                const satst_t *sa=bl->ss+cs[a].sat-1;
                int r=cs[a].sat;
                nk=rtk?o->nf:1;
                for (k=0;k<nk;k++) {
                    double N,Nr;
                    if (k==2&&!(st->ewl&&sa->ewl&&rtk_bias_ok(&bl->rtk,s,2)&&rtk_bias_ok(&bl->rtk,r,2))) continue;
                    N =k==0?st->N1:(k==1?N2of(st):N3of(st));
                    Nr=k==0?sa->N1:(k==1?N2of(sa):N3of(sa));
                    c[nc].n=2;
                    c[nc].i[0]=rtk?IB_(s,k,o):s; c[nc].h[0]= 1.0;
                    c[nc].i[1]=rtk?IB_(r,k,o):r; c[nc].h[1]=-1.0;
                    c[nc].v=(N-Nr)-(x[c[nc].i[0]]-x[c[nc].i[1]]);
                    c[nc].var=var_nl;
                    nc++;
                }
            }
            if (w>=0&&i!=w&&st->wl&&rtk_bias_ok(&bl->rtk,s,0)&&rtk_bias_ok(&bl->rtk,s,1)) {
                const satst_t *sw=bl->ss+cs[w].sat-1;
                int r=cs[w].sat;
                if (!(st->nl==1&&sw->nl==1)) {
                    c[nc].n=4;
                    c[nc].i[0]=IB_(s,0,o); c[nc].h[0]= 1.0;
                    c[nc].i[1]=IB_(s,1,o); c[nc].h[1]=-1.0;
                    c[nc].i[2]=IB_(r,0,o); c[nc].h[2]=-1.0;
                    c[nc].i[3]=IB_(r,1,o); c[nc].h[3]= 1.0;
                    c[nc].v=(st->Nwl-sw->Nwl)-((x[c[nc].i[0]]-x[c[nc].i[1]])-
                                              (x[c[nc].i[2]]-x[c[nc].i[3]]));
                    c[nc].var=var_wl;
                    nc++;
                }
            }
            if (e>=0&&i!=e&&st->ewl&&rtk_bias_ok(&bl->rtk,s,1)&&rtk_bias_ok(&bl->rtk,s,2)) {
                const satst_t *se=bl->ss+cs[e].sat-1;
                int r=cs[e].sat;
                if (!(st->nl==1&&se->nl==1&&a>=0)) {
                    c[nc].n=4;
                    c[nc].i[0]=IB_(s,1,o); c[nc].h[0]= 1.0;
                    c[nc].i[1]=IB_(s,2,o); c[nc].h[1]=-1.0;
                    c[nc].i[2]=IB_(r,1,o); c[nc].h[2]=-1.0;
                    c[nc].i[3]=IB_(r,2,o); c[nc].h[3]= 1.0;
                    c[nc].v=(st->Newl-se->Newl)-((x[c[nc].i[0]]-x[c[nc].i[1]])-
                                                (x[c[nc].i[2]]-x[c[nc].i[3]]));
                    c[nc].var=var_wl;
                    nc++;
                }
            }
        }
    }
    return nc;
}
static int apply_cons(int nx, double *x, double *P, const cons_t *c, int nc)
{
    double *H,*v,*R;
    int i,k,info;

    if (nc<=0) return 0;
    H=zeros(nx,nc); v=mat(nc,1); R=zeros(nc,nc);
    for (k=0;k<nc;k++) {
        for (i=0;i<c[k].n;i++) H[c[k].i[i]+k*nx]=c[k].h[i];
        v[k]=c[k].v;
        R[k+k*nc]=c[k].var;
    }
    info=filter(x,P,H,v,R,nx,nc);
    free(H); free(v); free(R);
    return info;
}
/* wide-lane / extra-wide-lane fixing from smoothed MW ------------------------
 * The SD MW mean of satellite s, referred to every already fixed satellite k
 * of the system, gives the float SD wide-lane  x_k = mw(s)-mw(k)+N(k). The
 * weighted mean over all k is rounded, so that multipath on a single
 * reference satellite cannot bias the result. Fixed wide-lanes are
 * re-checked the same way.                                                   */
static int mw_ready(const nrtk_baseline_t *bl, const satst_t *st, int k, gtime_t t)
{
    double tmin=k==0?bl->opt.wl_min_time:0.3*bl->opt.wl_min_time;
    double nmin=k==0?bl->opt.wl_min_epochs:MAX(3.0,0.5*bl->opt.wl_min_epochs);
    return st->nmw[k]>=nmin&&st->mwt0[k].time&&timediff(t,st->mwt0[k])>=tmin-1E-6;
}
static int wl_estimate(const nrtk_baseline_t *bl, const comsat_t *cs, int n, int sys,
                       int k, int self, double *x, double *sd)
{
    const satst_t *st=bl->ss+cs[self].sat-1;
    double sw=0.0,swx=0.0;
    int i,cnt=0;

    for (i=0;i<n;i++) {
        const satst_t *sk=bl->ss+cs[i].sat-1;
        double v,w;
        if (i==self||cs[i].sysidx!=sys||!(k==0?sk->wl:sk->ewl)) continue;
        if (!(k==0?cs[i].ok2:cs[i].ok3)||sk->nmw[k]<2) continue;
        v=mw_var_mean(st,k)+mw_var_mean(sk,k);
        w=1.0/v;
        swx+=w*(mw_mean(st,k)-mw_mean(sk,k)+(k==0?sk->Nwl:sk->Newl));
        sw+=w; cnt++;
    }
    if (cnt==0) return 0;
    *x=swx/sw;
    *sd=sqrt(1.0/sw);
    return cnt;
}
static void fix_widelanes(nrtk_baseline_t *bl, gtime_t t, const comsat_t *cs, int n)
{
    int sys,k,i,best,nfix;
    double x,sd;

    for (k=0;k<2;k++) {
        for (sys=0;sys<NRTK_NSYS;sys++) {
            /* re-check fixed wide-lanes against the other fixed ones */
            for (i=0;i<n;i++) {
                satst_t *st=bl->ss+cs[i].sat-1;
                if (cs[i].sysidx!=sys||!(k==0?st->wl:st->ewl)) continue;
                if (st->nmw[k]<30) continue;
                if (wl_estimate(bl,cs,n,sys,k,i,&x,&sd)<2) continue;
                if (fabs(x-(k==0?st->Nwl:st->Newl))>0.45) {
                    nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: %s x=%.2f N=%.0f",
                             bl->ista_m,bl->ista_a,cs[i].sat,k?"EWL":"WL",x,
                             k==0?st->Nwl:st->Newl);
                    if (k==0) {
                        /* the L1 integer was derived from this wide-lane and the
                           float is pinned to it: start that satellite over */
                        if (st->nl) {bl->nrefix++; reset_amb(bl,cs[i].sat);}
                        st->wl=0; st->nl=0;
                    }
                    else st->ewl=0;
                    reset_mw(st,k);
                }
            }
            /* datum: no fixed satellite in the system yet */
            for (i=nfix=0,best=-1;i<n;i++) {
                const satst_t *st=bl->ss+cs[i].sat-1;
                if (cs[i].sysidx!=sys||!(k==0?cs[i].ok2:cs[i].ok3)) continue;
                if (k==0?st->wl:st->ewl) {nfix++; continue;}
                if (mw_ready(bl,st,k,t)&&cs[i].el>=bl->opt.elmin_ar&&
                    (best<0||cs[i].el>cs[best].el)) best=i;
            }
            if (nfix==0) {
                if (best<0) continue;
                satst_t *st=bl->ss+cs[best].sat-1;   /* arbitrary integer */
                if (k==0) {st->Nwl=ROUND(mw_mean(st,0)); st->wl=1;}
                else      {st->Newl=ROUND(mw_mean(st,1)); st->ewl=1;}
            }
            for (i=0;i<n;i++) {
                satst_t *st=bl->ss+cs[i].sat-1;
                if (cs[i].sysidx!=sys||!(k==0?cs[i].ok2:cs[i].ok3)) continue;
                if ((k==0?st->wl:st->ewl)||!mw_ready(bl,st,k,t)) continue;
                if (!wl_estimate(bl,cs,n,sys,k,i,&x,&sd)) continue;
                if (fabs(x-ROUND(x))<=bl->opt.wl_max_frac&&sd<=bl->opt.wl_max_std*(k?2.0:1.0)) {
                    if (k==0) {st->Nwl=ROUND(x); st->wl=1;}
                    else      {st->Newl=ROUND(x); st->ewl=1;}
                }
            }
        }
    }
}
/* DD phase residual on proc slot k with SD integers -----------------------------*/
static double dd_phase_res(const comsat_t *s, const comsat_t *r, int k,
                           double Ns, double Nr)
{
    double lam=CLIGHT/s->m->freq[k],lamr=CLIGHT/r->m->freq[k];
    double sds=(nrtk_phase_m(s->a,k)-s->a->r)-(nrtk_phase_m(s->m,k)-s->m->r)-lam*Ns;
    double sdr=(nrtk_phase_m(r->a,k)-r->a->r)-(nrtk_phase_m(r->m,k)-r->m->r)-lamr*Nr;
    return sds-sdr;
}
static double dd_nondisp(const comsat_t *s, const comsat_t *r, const satst_t *ss,
                         const satst_t *sr, double *iono)
{
    double d1=dd_phase_res(s,r,0,ss->N1,sr->N1);
    double d2=dd_phase_res(s,r,1,N2of(ss),N2of(sr));
    double I,D;
    nrtk_split_iono(s->m->freq[0],s->m->freq[1],d1,d2,&I,&D);
    if (iono) *iono=I;
    return D;
}
/* troposphere model ----------------------------------------------------------------
 * SD (aux - master) Saastamoinen+GMF delay and station averaged wet mapping */
static void sd_trop(const nrtk_baseline_t *bl, gtime_t t, const comsat_t *c,
                    double *T, double *h)
{
    const double ep2000[]={2000,1,1,12,0,0};
    double pa[3],pm[3],mjd,mh,mwa,mwm;
    ecef2pos(bl->ra,pa); ecef2pos(bl->rm,pm);
    *T=nrtk_trop_slant(t,pa,c->a->azel)-nrtk_trop_slant(t,pm,c->m->azel);
    mjd=51544.5+timediff(t,epoch2time(ep2000))/86400.0;
    nrtk_gmf(mjd,pa[0],pa[1],pa[2],PI/2.0-c->a->azel[1],&mh,&mwa);
    nrtk_gmf(mjd,pm[0],pm[1],pm[2],PI/2.0-c->m->azel[1],&mh,&mwm);
    *h=0.5*(mwa+mwm);
}
/* independent model check: R = D - DD model must be explained by one relative
 * zenith wet delay  R ~ dzwd*h. A wrong L1 integer shifts D by lam_NL
 * (~0.11 m), a wrong wide-lane by lam_NL*f2/(f1-f2) (~0.32-0.38 m). */
static void model_res(const nrtk_baseline_t *bl, gtime_t t, const comsat_t *s,
                      const comsat_t *r, double D, double *R, double *h)
{
    double Ts,Tr,hs,hr;
    sd_trop(bl,t,s,&Ts,&hs);
    sd_trop(bl,t,r,&Tr,&hr);
    *R=D-(Ts-Tr);
    *h=hs-hr;
}
/* robust weighted LS of dzwd with prior, outliers beyond max(3 sigma, lim) */
static int fit_dzwd(const double *R, const double *h, const double *w, int n,
                    double lim, double *dz)
{
    double swh=0,swr=0,e[MAXOBS],med,mad,thr;
    int i,j,k,iter,use[MAXOBS],nu=0;

    if (n<=0) return 0;
    for (i=0;i<n;i++) use[i]=1;
    for (iter=0;iter<3;iter++) {
        swh=swr=0.0;
        for (i=nu=0;i<n;i++) if (use[i]) {swh+=w[i]*h[i]*h[i]; swr+=w[i]*h[i]*R[i]; nu++;}
        if (nu==0) return 0;
        /* prior dzwd ~ N(0,DZWD_PRIOR^2) relative to the weights (phase sigma 1 cm) */
        *dz=swr/(swh+SQR(0.01/DZWD_PRIOR));
        for (i=k=0;i<n;i++) if (use[i]) e[k++]=fabs(R[i]-*dz*h[i]);
        for (i=0;i<k;i++) for (j=i+1;j<k;j++) if (e[j]<e[i]) {double x=e[i]; e[i]=e[j]; e[j]=x;}
        med=e[k/2]; mad=1.4826*med;
        thr=MAX(3.0*mad,lim);
        for (i=0;i<n;i++) use[i]=fabs(R[i]-*dz*h[i])<=thr;
    }
    return nu;
}
static void update_dzwd(nrtk_baseline_t *bl, gtime_t t, const comsat_t *cs, int n)
{
    double R[MAXOBS],h[MAXOBS],w[MAXOBS],D,dz;
    int sys,i,m=0,a;

    for (sys=0;sys<NRTK_NSYS;sys++) {
        if ((a=anchor(bl,cs,n,sys,2))<0) continue;
        for (i=0;i<n&&m<MAXOBS;i++) {
            const satst_t *st=bl->ss+cs[i].sat-1;
            if (i==a||cs[i].sysidx!=sys||st->nl!=1||st->newfix||!cs[i].ok2) continue;
            D=dd_nondisp(cs+i,cs+a,st,bl->ss+cs[a].sat-1,NULL);
            model_res(bl,t,cs+i,cs+a,D,R+m,h+m);
            w[m]=SQR(sin(cs[i].el));
            m++;
        }
    }
    if (m>=3&&fit_dzwd(R,h,w,m,bl->opt.max_res,&dz)>=3&&fabs(dz)<=DZWD_MAX) {
        bl->dzwd=bl->dzwd_ok?0.8*bl->dzwd+0.2*dz:dz;
        bl->dzwd_ok=1;
    }
    if (bl->nep%60==0&&nrtk_log_enabled(NRTK_LOG_DEBUG)) {
        nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d dzwd=%.3f ok=%d (fit %.3f from %d sats)",
                 bl->ista_m,bl->ista_a,bl->dzwd,bl->dzwd_ok,m>=3?dz:0.0,m);
    }
}
static double wl_step(const comsat_t *s)   /* D shift per wide-lane cycle (m) */
{
    double f1=s->m->freq[0],f2=s->m->freq[1];
    return CLIGHT/(f1+f2)*f2/(f1-f2);
}
/* validate tentatively fixed satellites (nl==2) --------------------------------
 * A wrong integer of a partner shows up in its own residual only; a wrong
 * reference shifts every partner by the same amount. When all (>=2) partners
 * fail with a common offset the reference is suspect: a new (nl==2) datum is
 * dropped at once, a held one after two such rounds, and the partners are not
 * penalized (otherwise a wrong held reference would reject every new fix of
 * its system until it sets). */
static void validate_new(nrtk_baseline_t *bl, gtime_t t, const comsat_t *cs, int n,
                         const int *ref)
{
    double R[MAXOBS],h[MAXOBS],w[MAXOBS],e[MAXOBS],dz=0.0,k,emin,emax;
    int idx[MAXOBS],sys,i,m=0,bad,cnt,refbad;

    for (sys=0;sys<NRTK_NSYS;sys++) {
        int r=ref[sys];
        if (r<0||!bl->ss[cs[r].sat-1].nl||!cs[r].ok2) continue;
        for (i=0;i<n&&m<MAXOBS;i++) {
            satst_t *st=bl->ss+cs[i].sat-1;
            if (i==r||cs[i].sysidx!=sys||!st->nl||!cs[i].ok2) continue;
            model_res(bl,t,cs+i,cs+r,dd_nondisp(cs+i,cs+r,st,bl->ss+cs[r].sat-1,NULL),R+m,h+m);
            w[m]=SQR(sin(cs[i].el));
            idx[m++]=i;
        }
    }
    /* relative ZWD: from held satellites if known, else from this batch */
    if (bl->dzwd_ok) dz=bl->dzwd;
    else if (!fit_dzwd(R,h,w,m,bl->opt.max_res,&dz)) dz=0.0;
    if (fabs(dz)>DZWD_MAX) dz=0.0;   /* implausible: test against the plain model */

    for (sys=0;sys<NRTK_NSYS;sys++) {
        int r=ref[sys];
        satst_t *sr;
        if (r<0||!bl->ss[cs[r].sat-1].nl||!cs[r].ok2) continue;
        sr=bl->ss+cs[r].sat-1;
        emin=1E9; emax=-1E9;
        for (i=bad=cnt=0;i<m;i++) {
            int j=idx[i];
            if (cs[j].sysidx!=sys||bl->ss[cs[j].sat-1].nl!=2) continue;
            cnt++;
            e[i]=R[i]-dz*h[i];
            if (fabs(e[i])<=bl->opt.max_res) continue;
            bad++;
            emin=MIN(emin,e[i]); emax=MAX(emax,e[i]);
        }
        if (cnt==0) continue;
        /* every partner fails by (nearly) the same amount: the reference is suspect */
        refbad=cnt>=2&&bad==cnt&&emax-emin<=2.0*bl->opt.max_res;
        if (sr->nl==2) {
            if (bad==cnt) {            /* no partner accepted: drop the new datum */
                sr->nl=0;
                if (refbad) {sr->wl=0; reset_mw(sr,0);}
            }
        }
        else if (refbad) {
            if (++sr->nrej>=2) {
                nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: reference fails validation "
                         "(%d partners, res %.3f..%.3f)",bl->ista_m,bl->ista_a,cs[r].sat,cnt,emin,emax);
                sr->nl=0; sr->nrej=0; sr->wl=0;
                reset_mw(sr,0); reset_amb(bl,cs[r].sat);
                bl->nrefix++;
            }
        }
        else if (bad<cnt) sr->nrej=0;

        for (i=0;i<m;i++) {
            int j=idx[i];
            satst_t *st=bl->ss+cs[j].sat-1;
            if (cs[j].sysidx!=sys||st->nl!=2||fabs(e[i])<=bl->opt.max_res) continue;
            st->nl=0;
            if (refbad) continue;        /* not the partner's fault */
            k=ROUND(e[i]/wl_step(cs+j));
            nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d fix rejected res=%.3f dzwd=%.3f",
                     bl->ista_m,bl->ista_a,cs[j].sat,e[i],dz);
            if ((k!=0.0&&fabs(e[i]-k*wl_step(cs+j))<bl->opt.max_res)||++st->nrej>=3) {
                nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: WL suspect (res=%.3f k=%.0f)",
                         bl->ista_m,bl->ista_a,cs[j].sat,e[i],k);
                st->wl=0; st->nrej=0;
                reset_mw(st,0);
            }
        }
    }
}
/* monitor held satellites against the model ---------------------------------------
 * Per epoch: |e| > 2*max_res for 5 epochs (wrong wide-lane, slips). Averaged:
 * a wrong L1 integer shifts e by lam_NL (~0.107 m); fixed at low elevation where
 * multipath and troposphere residuals reach several cm it may pass the
 * validation, so the mean residual is tested again once the satellite is high. */
static void monitor_held(nrtk_baseline_t *bl, gtime_t t, const comsat_t *cs, int n)
{
    double R,h,e,dt;
    int sys,i,a,nbad,nfix,bad[MAXOBS];

    if (!bl->dzwd_ok) return;
    for (sys=0;sys<NRTK_NSYS;sys++) {
        const satst_t *sa;
        if ((a=anchor(bl,cs,n,sys,2))<0) continue;
        sa=bl->ss+cs[a].sat-1;
        for (i=nbad=nfix=0;i<n;i++) {
            satst_t *st=bl->ss+cs[i].sat-1;
            bad[i]=0;
            if (i==a||cs[i].sysidx!=sys||st->nl!=1||!cs[i].ok2) continue;
            nfix++;
            model_res(bl,t,cs+i,cs+a,dd_nondisp(cs+i,cs+a,st,bl->ss+cs[a].sat-1,NULL),&R,&h);
            e=R-bl->dzwd*h;
            if (fabs(e)>2.0*bl->opt.max_res) st->nbad++; else st->nbad=0;

            dt=st->eavg_t.time?timediff(t,st->eavg_t):1E9;
            if (st->eavg_w>0.0&&dt<=30.0&&st->eavg_anc!=cs[a].sat) {
                /* anchor changed: re-reference the mean with the new anchor's
                   own mean w.r.t. the old one, DD(i,new) = DD(i,old) - DD(new,old),
                   instead of restarting the average (a wrong old anchor cancels) */
                if (sa->eavg_w>0.0&&sa->eavg_anc==st->eavg_anc&&sa->eavg_t.time&&
                    timediff(t,sa->eavg_t)<=30.0) {
                    st->eavg-=sa->eavg;
                    st->eavg_w=MIN(st->eavg_w,sa->eavg_w);
                    st->eavg_anc=cs[a].sat;
                }
                else st->eavg_w=0.0;
            }
            if (st->eavg_w<=0.0||dt>30.0) {
                st->eavg=e; st->eavg_w=1E-9; st->eavg_anc=cs[a].sat;
            }
            else {
                st->eavg+=(1.0-exp(-dt/HELD_TAU))*(e-st->eavg);
                st->eavg_w+=dt;
            }
            st->eavg_t=t;
            bad[i]=st->nbad>=5?1:
                   (cs[i].el>=HELD_EL&&st->eavg_w>=HELD_MIN&&fabs(st->eavg)>HELD_THR?2:0);
            if (bad[i]) nbad++;
        }
        if (nfix>=2&&nbad*2>nfix) {          /* anchor suspect */
            satst_t *sa=bl->ss+cs[a].sat-1;
            nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: anchor fails model",
                     bl->ista_m,bl->ista_a,cs[a].sat);
            sa->nl=0; sa->wl=0; reset_mw(sa,0); reset_amb(bl,cs[a].sat); bl->nrefix++;
            for (i=0;i<n;i++) {
                if (cs[i].sysidx!=sys) continue;
                bl->ss[cs[i].sat-1].nbad=0;
                bl->ss[cs[i].sat-1].eavg_w=0.0;
            }
            continue;
        }
        for (i=0;i<n;i++) {
            satst_t *st=bl->ss+cs[i].sat-1;
            if (i==a||cs[i].sysidx!=sys||st->nl!=1||!bad[i]) continue;
            if (bad[i]==1) {
                nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: fails model",
                         bl->ista_m,bl->ista_a,cs[i].sat);
                st->wl=0; reset_mw(st,0);
            }
            else {
                nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: mean model residual %.3f m",
                         bl->ista_m,bl->ista_a,cs[i].sat,st->eavg);
            }
            st->nl=0; st->nbad=0; reset_amb(bl,cs[i].sat); bl->nrefix++;
        }
    }
}
/* rtklib engine -------------------------------------------------------------------*/
static int build_obs(const nrtk_baseline_t *bl, const nrtk_staepoch_t *em,
                     const nrtk_staepoch_t *ea, const comsat_t *cs,
                     const int *slip, int n, obsd_t *obs)
{
    int i,j,k,nobs=0,r;
    const nrtk_staepoch_t *ep[2]={ea,em};

    for (r=0;r<2;r++) {
        for (i=0;i<ep[r]->ns;i++) {
            const nrtk_satdata_t *s=ep[r]->sat+i;
            obsd_t *o;
            if (!(bl->opt.sysmask&(1<<s->sysidx))) continue;
            if (nobs>=MAXOBS*2) break;
            o=obs+nobs++;
            *o=s->obs;
            o->rcv=(uint8_t)(r+1);
            o->time=ep[0]->time;
            for (j=bl->nf;j<NFREQ+NEXOBS;j++) {
                o->L[j]=o->P[j]=0.0; o->D[j]=0.0f; o->code[j]=CODE_NONE; o->LLI[j]=0;
            }
            for (k=0;k<n;k++) {
                if (cs[k].sat!=s->sat) continue;
                if (slip[k]) for (j=0;j<bl->nf;j++) o->LLI[j]|=LLI_SLIP;
                break;
            }
        }
    }
    return nobs;
}
static void model_zwd(const double *rr, double *zwd)
{
    double pos[3],zhd;
    ecef2pos(rr,pos);
    nrtk_trop_zenith(pos,0.7,&zhd,zwd);
}
static void set_state(double *x, double *P, int nx, int i, double xi, double var)
{
    int j;
    x[i]=xi;
    for (j=0;j<nx;j++) P[i+j*nx]=P[j+i*nx]=0.0;
    P[i*(nx+1)]=var;
}
/* Over short baselines only the relative zenith delay is observable in double
 * differences; a free common mode drifts and degrades low-elevation DDs. The
 * master (base) ZWD is pinned to the model every epoch and only the aux
 * (rover) ZWD is estimated, starting from its model value. */
static void rtk_trop(nrtk_baseline_t *bl)
{
    rtk_t *rtk=&bl->rtk;
    const prcopt_t *o=&rtk->opt;
    double zwd;
    int r,i;

    if (o->tropopt<TROPOPT_EST) return;
    for (r=bl->trop_init?1:0;r<2;r++) {     /* r: 0=aux(rover), 1=master(base) */
        model_zwd(r?bl->rm:bl->ra,&zwd);
        set_state(rtk->x,rtk->P,rtk->nx,IT_(r,o),zwd,r?SQR(1E-3):SQR(REL_ZWD));
        if (o->tropopt>=TROPOPT_ESTG) {
            for (i=1;i<=2;i++) set_state(rtk->x,rtk->P,rtk->nx,IT_(r,o)+i,1E-6,r?SQR(1E-5):SQR(0.001));
        }
    }
    bl->trop_init=1;
}
static int rtk_update(nrtk_baseline_t *bl, const nrtk_staepoch_t *em,
                      const nrtk_staepoch_t *ea, const comsat_t *cs,
                      const int *slip, int n, const nav_t *nav)
{
    obsd_t *obs=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS*2);
    int nobs,stat;

    rtk_trop(bl);
    nobs=build_obs(bl,em,ea,cs,slip,n,obs);
    stat=rtkpos(&bl->rtk,obs,nobs,nav)&&bl->rtk.sol.stat!=SOLQ_NONE;
    free(obs);
    return stat;
}
/* native engine -------------------------------------------------------------------
 * SD ionosphere-free observable minus geometry, model troposphere and the
 * fixed wide-lane term:  y = SD[Phi_IF - rho - T] - c*f2/(f1^2-f2^2)*Nwl
 *                          = lam_NL*b + h*dzwd + (receiver terms) + e
 * DD within each system removes the receiver terms; b is the float L1
 * ambiguity (cycles) with the per-system datum of the SD integers. */
static double nat_obs(const nrtk_baseline_t *bl, gtime_t t, const comsat_t *c,
                      double Nwl, double *lam, double *h, double *var)
{
    double f1=c->m->freq[0],f2=c->m->freq[1],g,k,ifa,ifm,T,sinel=sin(c->el);
    g=f1*f1/(f1*f1-f2*f2); k=-f2*f2/(f1*f1-f2*f2);
    ifa=g*nrtk_phase_m(c->a,0)+k*nrtk_phase_m(c->a,1)-c->a->r;
    ifm=g*nrtk_phase_m(c->m,0)+k*nrtk_phase_m(c->m,1)-c->m->r;
    sd_trop(bl,t,c,&T,h);
    *lam=CLIGHT/(f1+f2);
    *var=(g*g+k*k)*2.0*(SQR(bl->opt.err_phase)+SQR(bl->opt.err_phase/sinel));
    return ifa-ifm-T-CLIGHT*f2/(f1*f1-f2*f2)*Nwl;
}
static void nat_clear(natf_t *f, int i)
{
    int j;
    f->x[i]=0.0;
    for (j=0;j<NAT_NX;j++) f->P[i+j*NAT_NX]=f->P[j+i*NAT_NX]=0.0;
}
static int nat_update(nrtk_baseline_t *bl, gtime_t t, const comsat_t *cs, int n,
                      const int *slip)
{
    natf_t *f=&bl->nat;
    const int nx=NAT_NX;
    double y[MAXOBS],lam[MAXOBS],h[MAXOBS],var[MAXOBS],*H,*v,*R,dt;
    int use[MAXOBS],rej[MAXOBS],row_i[MAXOBS],row_p[MAXOBS],i,j,k,s,sys,m=0,p,ps;

    dt=f->t.time?timediff(t,f->t):0.0;
    if (!f->init) {set_state(f->x,f->P,nx,0,1E-6,SQR(REL_ZWD)); f->init=1;}
    else if (dt>0.0) f->P[0]+=SQR(bl->opt.trop_prn)*dt;

    /* usable satellites: wide-lane fixed, both frequencies, no slip */
    for (i=0;i<n;i++) {
        satst_t *st=bl->ss+cs[i].sat-1;
        use[i]=st->wl&&cs[i].ok2&&!slip[i];
        if (!use[i]) continue;
        y[i]=nat_obs(bl,t,cs+i,st->Nwl,lam+i,h+i,var+i);
    }
    /* drop the states of satellites that slipped or whose wide-lane changed; a
       satellite merely missing an observation this epoch (or temporarily below
       the mask) keeps its state and is just not updated (gaps > MAXGAP are
       slips, see detect_slips) */
    for (s=1;s<=MAXSAT;s++) {
        if (f->x[s]==0.0) continue;
        for (i=0;i<n;i++) if (cs[i].sat==s) break;
        if (i<n&&(!bl->ss[s-1].wl||slip[i]||bl->ss[s-1].Nwl_nat!=bl->ss[s-1].Nwl)) nat_clear(f,s);
        else f->P[s*(nx+1)]+=SQR(1E-4)*fabs(dt);
    }
    H=zeros(nx,n); v=mat(n,1);
    for (sys=0;sys<NRTK_NSYS;sys++) {
        double pred,vv,hh,q;
        int nbad=0,ntest=0,m0=m;
        /* pivot: highest elevation satellite with an active state */
        for (i=0,p=-1;i<n;i++) {
            if (!use[i]||cs[i].sysidx!=sys) continue;
            if (f->x[cs[i].sat]!=0.0&&(p<0||cs[i].el>cs[p].el)) p=i;
        }
        if (p<0) {
            for (i=0;i<n;i++) {
                if (!use[i]||cs[i].sysidx!=sys) continue;
                if (p<0||cs[i].el>cs[p].el) p=i;
            }
            if (p<0) continue;
            ps=cs[p].sat;
            set_state(f->x,f->P,nx,ps,y[p]/lam[p],VAR_NEWAMB);
            if (f->x[ps]==0.0) f->x[ps]=1E-9;
            bl->ss[ps-1].Nwl_nat=bl->ss[ps-1].Nwl;
        }
        ps=cs[p].sat;
        for (i=0;i<n;i++) {
            s=cs[i].sat;
            rej[i]=0;
            if (!use[i]||cs[i].sysidx!=sys||i==p) continue;
            if (f->x[s]==0.0) {   /* new state consistent with the DD measurement */
                set_state(f->x,f->P,nx,s,(y[i]-y[p]+lam[p]*f->x[ps]-(h[i]-h[p])*f->x[0])/lam[i],
                          VAR_NEWAMB);
                if (f->x[s]==0.0) f->x[s]=1E-9;
                bl->ss[s-1].Nwl_nat=bl->ss[s-1].Nwl;
            }
            hh=h[i]-h[p];
            pred=lam[i]*f->x[s]-lam[p]*f->x[ps]+hh*f->x[0];
            vv=(y[i]-y[p])-pred;
            q=SQR(lam[i])*f->P[s*(nx+1)]+SQR(lam[p])*f->P[ps*(nx+1)]-
              2.0*lam[i]*lam[p]*f->P[s+ps*nx]+SQR(hh)*f->P[0]+
              2.0*hh*(lam[i]*f->P[s*nx]-lam[p]*f->P[ps*nx]);
            if (q<SQR(0.05)) ntest++;          /* converged: the test is meaningful */
            if (fabs(vv)>MAX(0.25,5.0*sqrt(q+var[i]+var[p]))) {
                rej[i]=1; nbad++;
                continue;
            }
            H[0+m*nx]=hh; H[s+m*nx]=lam[i]; H[ps+m*nx]=-lam[p];
            v[m]=vv; row_i[m]=i; row_p[m]=p;
            m++;
        }
        /* an outlier is an undetected slip: of the satellite itself, or of the
           pivot when most converged partners fail together. Either way the
           satellite starts over (the wide-lane may have changed as well) */
        if (nbad>=2&&2*nbad>ntest) {
            nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d native pivot outlier (%d of %d partners)",
                     bl->ista_m,bl->ista_a,ps,nbad,ntest);
            if (bl->ss[ps-1].nl) bl->nrefix++;
            reset_sat(bl->ss+ps-1); nat_clear(f,ps); bl->nslip++;
            memset(H+m0*nx,0,sizeof(double)*nx*(m-m0));   /* this epoch's rows of the system */
            m=m0;
            continue;
        }
        for (i=0;i<n;i++) {
            if (!rej[i]) continue;
            s=cs[i].sat;
            nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d native outlier",bl->ista_m,bl->ista_a,s);
            if (bl->ss[s-1].nl) bl->nrefix++;
            reset_sat(bl->ss+s-1); nat_clear(f,s); bl->nslip++;
        }
    }
    if (m>0) {
        R=zeros(m,m);
        for (j=0;j<m;j++) for (k=0;k<m;k++) {
            if (j==k) R[j+k*m]=var[row_i[j]]+var[row_p[j]];
            else if (row_p[j]==row_p[k]) R[j+k*m]=var[row_p[j]];
        }
        if (f->x[0]==0.0) f->x[0]=1E-12;
        if (filter(f->x,f->P,H,v,R,nx,m)) {
            nrtk_log(NRTK_LOG_WARN,"netpos","bl %d-%d native filter error",bl->ista_m,bl->ista_a);
        }
        free(R);
    }
    free(H); free(v);
    f->t=t;
    return m>0;
}
/* L1 ambiguity resolution ---------------------------------------------------------*/
static int resolve_l1(nrtk_baseline_t *bl, gtime_t t, const comsat_t *cs, int n)
{
    int nx=eng_nx(bl),i,j,k,m,nc,nfixed=0,ncand=0,cand[MAXOBS],ref[NRTK_NSYS];
    int idx[MAXOBS],rix[MAXOBS],cidx[MAXOBS],drop,info,nb,accepted=0;
    double *xc,*Pc,*y,*Q,*F,s[2];
    cons_t *c;

    for (i=0;i<n;i++) if (bl->ss[cs[i].sat-1].nl) nfixed++;

    /* candidates: wide-lane fixed, locked, above AR mask, states valid */
    for (i=0;i<n;i++) {
        const satst_t *st=bl->ss+cs[i].sat-1;
        if (st->nl||!st->wl||!cs[i].ok2) continue;
        if (st->lock<bl->opt.min_lock||cs[i].el<bl->opt.elmin_ar) continue;
        if (!amb_ok(bl,cs[i].sat)) continue;
        if (bl->engine==NRTK_ENGINE_RTKLIB&&
            (!bl->rtk.ssat[cs[i].sat-1].vsat[0]||!bl->rtk.ssat[cs[i].sat-1].vsat[1])) continue;
        cand[ncand++]=i;
    }
    if (ncand==0) return 0;
    for (i=0;i<ncand;i++) for (j=i+1;j<ncand;j++) {   /* elevation descending */
        if (cs[cand[j]].el>cs[cand[i]].el) {k=cand[i]; cand[i]=cand[j]; cand[j]=k;}
    }
    /* condition on the held integers (and, rtklib, the fixed wide-lanes) */
    xc=mat(nx,1); Pc=mat(nx,nx);
    matcpy(xc,eng_x(bl),nx,1); matcpy(Pc,eng_P(bl),nx,nx);
    c=(cons_t *)malloc(sizeof(cons_t)*MAXOBS*4);
    nc=collect_cons(bl,cs,n,VAR_COND,VAR_COND,1,c,MAXOBS*4);
    if (apply_cons(nx,xc,Pc,c,nc)) {free(xc); free(Pc); free(c); return 0;}
    free(c);

    y=mat(MAXOBS,1); Q=mat(MAXOBS,MAXOBS); F=mat(MAXOBS,2);

    for (drop=0;drop<=bl->opt.max_drop&&ncand>0;drop++) {
        /* DD reference per system: held anchor, else best candidate */
        for (m=0;m<NRTK_NSYS;m++) {
            ref[m]=anchor(bl,cs,n,m,2);
            if (ref[m]<0) {
                for (k=0;k<ncand;k++) if (cs[cand[k]].sysidx==m) {ref[m]=cand[k]; break;}
            }
        }
        for (k=nb=0;k<ncand;k++) {
            int ci=cand[k],r=ref[cs[ci].sysidx];
            double q;
            if (r<0||r==ci) continue;
            idx[nb]=amb_idx(bl,cs[ci].sat); rix[nb]=amb_idx(bl,cs[r].sat);
            q=Pc[idx[nb]*(nx+1)]-2.0*Pc[idx[nb]+rix[nb]*nx]+Pc[rix[nb]*(nx+1)];
            if (q>SQR(NL_MAXSTD)) continue;   /* float not converged */
            cidx[nb]=k;
            nb++;
        }
        if (nb<1||(nfixed==0&&nb<bl->opt.min_fix)) break;
        for (i=0;i<nb;i++) {
            y[i]=xc[idx[i]]-xc[rix[i]];
            for (j=0;j<nb;j++) {
                Q[i+j*nb]=Pc[idx[i]+idx[j]*nx]-Pc[idx[i]+rix[j]*nx]-
                          Pc[rix[i]+idx[j]*nx]+Pc[rix[i]+rix[j]*nx];
            }
        }
        info=lambda(nb,2,y,Q,F,s);
        bl->ratio=(!info&&s[0]>0.0)?MIN(s[1]/s[0],999.9):0.0;
        if (nrtk_log_enabled(NRTK_LOG_DEBUG)) {
            char buf[1024]; int p=0;
            for (i=0;i<nb&&p<900;i++) {
                p+=snprintf(buf+p,sizeof(buf)-p," %d:%.2f(%.2f)",amb_sat(bl,idx[i]),y[i],
                            sqrt(Q[i+i*nb]));
            }
            nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d AR cand=%d nb=%d ratio=%.2f%s",bl->ista_m,
                     bl->ista_a,ncand,nb,bl->ratio,buf);
        }
        if (!info&&(s[0]<=0.0||s[1]/s[0]>=bl->opt.ar_ratio)) {
            /* accept: SD integers (datum for systems without held satellites) */
            for (m=0;m<NRTK_NSYS;m++) {
                int r=ref[m];
                satst_t *sr;
                if (r<0) continue;
                sr=bl->ss+cs[r].sat-1;
                if (!sr->nl) {
                    int used=0;
                    for (i=0;i<nb;i++) if (rix[i]==amb_idx(bl,cs[r].sat)) used=1;
                    if (!used) continue;
                    sr->N1=ROUND(xc[amb_idx(bl,cs[r].sat)]);
                    sr->nl=2;   /* tentative */
                }
            }
            for (i=0;i<nb;i++) {
                satst_t *st=bl->ss+amb_sat(bl,idx[i])-1,*sr=bl->ss+amb_sat(bl,rix[i])-1;
                st->N1=sr->N1+F[i];
                st->nl=2;
            }
            validate_new(bl,t,cs,n,ref);
            for (i=0;i<n;i++) {
                satst_t *st=bl->ss+cs[i].sat-1;
                if (st->nl==2) {st->nl=1; st->newfix=1; st->nrej=0; accepted++;}
            }
            break;
        }
        /* partial fixing: remove one satellite of the search set. (Dropping the
           lowest candidate overall, as before, could remove satellites that were
           not even in the set because their float had not converged, leaving the
           set and the ratio unchanged.) Which one depends on the engine:
             native : the lowest elevation satellite. Its ionosphere-free float
                      is tightly converged, so a non-integer fraction reflects a
                      bias (multipath, troposphere) rather than noise, and the
                      lowest satellite carries the largest one. Dropping the DD
                      farthest from its integer instead fixed more satellites at
                      10-13 deg and raised the wrong-fix rate of the HKSC rover
                      from 0.12% to 1.3% (docs/VALIDATION.md).
             rtklib : a DD clearly off its integer (> PAR_AMBIG cycles) first -
                      with per-satellite ionosphere states the float is noisy
                      and such a DD blocks the ratio test whatever its elevation
                      (two Hong Kong baselines stayed at 0% and 31% otherwise);
                      the lowest elevation satellite when none is. */
        {
            int worst=-1,low=-1;
            double dmax=-1.0;
            for (i=0;i<nb;i++) {
                double d=info?0.0:fabs(y[i]-F[i]);
                if (d>dmax) {dmax=d; worst=i;}
                if (low<0||cidx[i]>cidx[low]) low=i;
            }
            if (low<0) break;
            if (bl->engine!=NRTK_ENGINE_RTKLIB||dmax<PAR_AMBIG) worst=low;
            for (k=cidx[worst];k<ncand-1;k++) cand[k]=cand[k+1];
            ncand--;
        }
    }
    free(y); free(Q); free(F); free(xc); free(Pc);
    return accepted;
}
/* hold: apply the held integers to the float filter --------------------------
 * Wide-lanes from MW only condition the search: a wrong wide-lane held in the
 * float filter would bias every later float solution. */
static void hold(nrtk_baseline_t *bl, const comsat_t *cs, int n)
{
    cons_t *c=(cons_t *)malloc(sizeof(cons_t)*MAXOBS*4);
    int nc=collect_cons(bl,cs,n,0.0,bl->opt.hold_var,0,c,MAXOBS*4);
    if (apply_cons(eng_nx(bl),eng_x(bl),eng_P(bl),c,nc)) {
        nrtk_log(NRTK_LOG_WARN,"netpos","bl %d-%d hold filter error",bl->ista_m,bl->ista_a);
    }
    free(c);
}
/* check held ambiguities against the float filter -----------------------------*/
static void check_held(nrtk_baseline_t *bl, const comsat_t *cs, int n)
{
    const double *x=eng_x(bl);
    int sys,i,a,k,nk,nbad,nfix,rtk=bl->engine==NRTK_ENGINE_RTKLIB;

    for (i=0;i<n;i++) {           /* float state reset */
        satst_t *st=bl->ss+cs[i].sat-1;
        if (st->nl&&!amb_ok(bl,cs[i].sat)) {
            nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: float state reset",
                     bl->ista_m,bl->ista_a,cs[i].sat);
            st->nl=0; bl->nrefix++;
        }
    }
    nk=rtk?2:1;
    for (sys=0;sys<NRTK_NSYS;sys++) {
        if ((a=anchor(bl,cs,n,sys,2))<0) continue;
        const satst_t *sa=bl->ss+cs[a].sat-1;
        for (i=nbad=nfix=0;i<n;i++) {
            satst_t *st=bl->ss+cs[i].sat-1;
            double dd;
            if (i==a||cs[i].sysidx!=sys||!st->nl) continue;
            nfix++;
            for (k=0;k<nk;k++) {
                dd=rtk?x[IB_(cs[i].sat,k,&bl->rtk.opt)]-x[IB_(cs[a].sat,k,&bl->rtk.opt)]:
                       x[cs[i].sat]-x[cs[a].sat];
                /* the float is closer to another integer. (0.5, not less: the
                   ratio test accepts floats up to ~0.37 cycles off, a tighter
                   limit here made such fixes flap every epoch) */
                if (fabs(dd-(k==0?st->N1-sa->N1:N2of(st)-N2of(sa)))>0.5) {
                    nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: held f%d dd=%.3f N=%.0f",
                             bl->ista_m,bl->ista_a,cs[i].sat,k,dd,k==0?st->N1-sa->N1:N2of(st)-N2of(sa));
                    st->nl=0; nbad++; bl->nrefix++;
                    break;
                }
            }
        }
        if (nfix>=2&&nbad*2>nfix) {  /* most partners inconsistent: anchor suspect */
            nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: anchor inconsistent",
                     bl->ista_m,bl->ista_a,cs[a].sat);
            bl->ss[cs[a].sat-1].nl=0;
        }
    }
}
/* cycle slip repair ------------------------------------------------------------------
 * A jump of the DD residuals of a held satellite by an integer number of cycles
 * on every frequency is a slip the detectors missed. Instead of starting the
 * satellite over, its SD integers, MW averages and float states are shifted by
 * the jump and the residuals stay continuous. */
static int jump_integers(const comsat_t *c, const double *dj, int has3, double *dn)
{
    int k,n=has3?3:2;
    for (k=0;k<n;k++) {
        double x=dj[k]/(CLIGHT/c->m->freq[k]);
        dn[k]=ROUND(x);
        if (fabs(x-dn[k])>0.15) return 0;
    }
    if (!has3) dn[2]=0.0;
    return dn[0]!=0.0||dn[1]!=0.0;
}
static void repair_slip(nrtk_baseline_t *bl, int sat, const double *dn, int has3)
{
    satst_t *st=bl->ss+sat-1;
    double dwl=dn[0]-dn[1],dewl=dn[1]-dn[2];
    int f,j;

    st->N1+=dn[0];
    st->Nwl+=dwl;
    if (dwl!=0.0) {         /* shift the MW sums: sum (x+d)^2 = sum x^2 + 2d sum x + n d^2 */
        st->mw2[0]+=2.0*dwl*st->mw[0]+st->nmw[0]*dwl*dwl;
        st->mw[0]+=st->nmw[0]*dwl;
    }
    if (st->ewl) {
        if (has3) {
            st->Newl+=dewl;
            if (dewl!=0.0) {
                st->mw2[1]+=2.0*dewl*st->mw[1]+st->nmw[1]*dewl*dewl;
                st->mw[1]+=st->nmw[1]*dewl;
            }
        }
        else {st->ewl=0; reset_mw(st,1);}  /* third frequency unknown */
    }
    st->Nwl_nat=st->Nwl;
    if (bl->engine==NRTK_ENGINE_RTKLIB) {
        for (f=0;f<bl->rtk.opt.nf&&f<3;f++) {
            if (f==2&&!has3) break;
            j=IB_(sat,f,&bl->rtk.opt);
            if (bl->rtk.x[j]!=0.0) bl->rtk.x[j]+=dn[f];
        }
    }
    else if (bl->nat.x[sat]!=0.0) bl->nat.x[sat]+=dn[0];   /* y shifts by lam_NL*dn1 */
}
/* statistics ---------------------------------------------------------------------*/
static void update_stat(nrtk_baseline_t *bl, gtime_t t, int fixed)
{
    if (bl->nep==0) bl->tfirst=t;
    bl->nep++;
    if (fixed) {
        bl->nfixep++;
        if (!bl->everfixed) {bl->everfixed=1; bl->ttff=timediff(t,bl->tfirst);}
    }
    bl->hist[bl->hi]=(uint8_t)fixed;
    bl->hi=(bl->hi+1)%60;
    if (bl->hn<60) bl->hn++;
    bl->fixed_prev=fixed;
}
/* process one epoch ----------------------------------------------------------------*/
extern int nrtk_bl_process(nrtk_baseline_t *bl, const nrtk_staepoch_t *em,
                           const nrtk_staepoch_t *ea, const nav_t *nav,
                           nrtk_blepoch_t *out)
{
    comsat_t cs[MAXOBS];
    double dj[MAXOBS][3];
    int slip[MAXOBS],jump[MAXOBS],tested[MAXOBS],has3[MAXOBS],oidx[MAXOBS];
    int i,k,n,sys,nfix=0,nwl=0,nslip=0,ok,ok3;

    memset(out,0,sizeof(*out));
    out->time=em->time;
    out->sta_m=bl->ista_m; out->sta_a=bl->ista_a;
    out->length=bl->length;
    out->dt_age=timediff(em->time,ea->time);

    if (fabs(out->dt_age)>DTTOL||em->ns<4||ea->ns<4) {
        update_stat(bl,em->time,0);
        return 0;
    }
    n=common_sats(bl,em,ea,cs);
    if (n<4) {update_stat(bl,em->time,0); return 0;}

    for (i=0;i<MAXSAT;i++) bl->ss[i].newfix=0;
    detect_slips(bl,em->time,cs,n,slip);
    for (i=0;i<n;i++) nslip+=slip[i];

    if (bl->engine==NRTK_ENGINE_RTKLIB) ok=rtk_update(bl,em,ea,cs,slip,n,nav);
    else {nat_update(bl,em->time,cs,n,slip); ok=1;}
    if (!ok) {update_stat(bl,em->time,0); return 0;}

    check_held(bl,cs,n);
    fix_widelanes(bl,em->time,cs,n);
    update_dzwd(bl,em->time,cs,n);
    monitor_held(bl,em->time,cs,n);
    resolve_l1(bl,em->time,cs,n);
    hold(bl,cs,n);

    /* outputs --------------------------------------------------------------------*/
    out->valid=1;
    out->ratio=bl->ratio;
    for (i=0;i<n;i++) {
        if (bl->ss[cs[i].sat-1].nl) nfix++;
        if (bl->ss[cs[i].sat-1].wl) nwl++;
    }
    for (sys=0;sys<NRTK_NSYS;sys++) {
        int r=anchor(bl,cs,n,sys,2);
        out->refsat[sys]=r>=0?cs[r].sat:0;
    }
    for (i=0;i<n;i++) {oidx[i]=-1; jump[i]=tested[i]=0;}
    for (i=0;i<n&&out->ns<NRTK_BL_MAXSAT;i++) {
        satst_t *st=bl->ss+cs[i].sat-1;
        nrtk_blsat_t *o=out->sat+out->ns;
        int r=-1;
        oidx[i]=out->ns++;
        o->sat=(uint8_t)cs[i].sat; o->sysidx=(int8_t)cs[i].sysidx;
        o->wlfixed=(uint8_t)st->wl; o->ewlfixed=(uint8_t)st->ewl;
        o->slip=(uint8_t)slip[i];
        o->azel_m[0]=(float)cs[i].m->azel[0]; o->azel_m[1]=(float)cs[i].m->azel[1];
        o->azel_a[0]=(float)cs[i].a->azel[0]; o->azel_a[1]=(float)cs[i].a->azel[1];
        if (st->nl) st->nfixep++; else st->nfixep=0;
        o->nfixep=st->nfixep;
        for (k=0;k<n;k++) if (cs[k].sat==out->refsat[cs[i].sysidx]) {r=k; break;}
        /* a held satellite without both phases this epoch has no residual to
           offer (it stays held, the float state is kept) */
        if (!st->nl||r<0||!cs[i].ok2) continue;
        o->fixed=1;
        o->newfix=(uint8_t)st->newfix;
        o->isref=(uint8_t)(r==i);
        if (r==i) continue;
        {
            const satst_t *sr=bl->ss+cs[r].sat-1;
            double Ns[3]={st->N1,N2of(st),N3of(st)},Nr[3]={sr->N1,N2of(sr),N3of(sr)};
            double pa[3],pm[3];
            for (k=0;k<bl->nf&&k<NRTK_NPROC;k++) {
                if (k==2&&!(st->ewl&&sr->ewl&&cs[i].ok3&&cs[r].ok3)) continue;
                o->dd_res[k]=dd_phase_res(cs+i,cs+r,k,Ns[k],Nr[k]);
                if (cs[i].a->obs.P[k]!=0.0&&cs[i].m->obs.P[k]!=0.0&&
                    cs[r].a->obs.P[k]!=0.0&&cs[r].m->obs.P[k]!=0.0) {
                    o->dd_code[k]=((cs[i].a->obs.P[k]-cs[i].a->r)-(cs[i].m->obs.P[k]-cs[i].m->r))-
                                  ((cs[r].a->obs.P[k]-cs[r].a->r)-(cs[r].m->obs.P[k]-cs[r].m->r));
                }
            }
            nrtk_split_iono(cs[i].m->freq[0],cs[i].m->freq[1],o->dd_res[0],o->dd_res[1],
                            &o->dd_iono,&o->dd_nondisp);
            ecef2pos(ea->rr,pa); ecef2pos(em->rr,pm);
            o->dd_trop_model=
                (nrtk_trop_slant(em->time,pa,cs[i].a->azel)-nrtk_trop_slant(em->time,pm,cs[i].m->azel))-
                (nrtk_trop_slant(em->time,pa,cs[r].a->azel)-nrtk_trop_slant(em->time,pm,cs[r].m->azel));
            /* continuity of the DD residuals of held satellites: an undetected
               slip or a wrong integer moves them by a wavelength */
            ok3=bl->nf>=3&&st->ewl&&sr->ewl&&cs[i].ok3&&cs[r].ok3;   /* dd_res[2] computed */
            if (st->res_ref==cs[r].sat&&!st->newfix) {
                tested[i]=1;
                has3[i]=ok3&&st->res_has3;
                for (k=0;k<3;k++) dj[i][k]=k<2||has3[i]?o->dd_res[k]-st->res_prev[k]:0.0;
                if (fabs(dj[i][0])>JUMP_RES||fabs(dj[i][1])>JUMP_RES) {
                    jump[i]=1;
                    nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d jump d1=%.3f d2=%.3f ref=%d",
                             bl->ista_m,bl->ista_a,cs[i].sat,dj[i][0],dj[i][1],cs[r].sat);
                }
            }
            for (k=0;k<3;k++) st->res_prev[k]=o->dd_res[k];
            st->res_has3=ok3;
            st->res_ref=cs[r].sat;
        }
    }
    /* a jump is a slip that escaped the detectors: of the satellite, or of the
       reference when most tested partners jump together. A jump by whole
       cycles on every frequency is repaired (integers, MW averages and float
       shifted; the residuals stay continuous). Anything else starts the
       satellite over like a detected slip: the wide-lane may have changed and
       the float is pinned to the old integer; for a reference the system's
       residuals of this epoch are unusable. */
    for (sys=0;sys<NRTK_NSYS;sys++) {
        int r=-1,nt=0,nj=0,h3=1;
        double dn[3],dnr[3];
        for (k=0;k<n;k++) if (cs[k].sat==out->refsat[sys]) {r=k; break;}
        if (r<0) continue;
        for (i=0;i<n;i++) if (cs[i].sysidx==sys&&tested[i]) {nt++; nj+=jump[i];}
        if (nj==0) continue;
        if (nj>=2&&2*nj>nt) {
            /* reference slipped: every tested partner jumped by the same integers
               (the third frequency from the partners that have it) */
            int agree=nj==nt,first=1,n3=0;
            dnr[0]=dnr[1]=dnr[2]=0.0;
            for (i=0;i<n&&agree;i++) {
                double neg[3];
                if (cs[i].sysidx!=sys||!jump[i]) continue;
                for (k=0;k<3;k++) neg[k]=-dj[i][k];
                if (!jump_integers(cs+i,neg,has3[i],dn)) {agree=0; break;}
                if (first) {dnr[0]=dn[0]; dnr[1]=dn[1]; first=0;}
                else if (dn[0]!=dnr[0]||dn[1]!=dnr[1]) agree=0;
                if (has3[i]) {
                    if (n3==0) dnr[2]=dn[2]; else if (dn[2]!=dnr[2]) agree=0;
                    n3++;
                }
            }
            h3=n3>0;
            if (agree) {
                nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d reference slip repaired (%.0f,%.0f,%.0f) from %d partners",
                         bl->ista_m,bl->ista_a,cs[r].sat,dnr[0],dnr[1],dnr[2],nj);
                repair_slip(bl,cs[r].sat,dnr,h3);
                bl->nslip++;
                if (oidx[r]>=0) {           /* VRS lock of the reference restarts */
                    out->sat[oidx[r]].slip=1;
                    out->sat[oidx[r]].ewlfixed=(uint8_t)bl->ss[cs[r].sat-1].ewl;
                }
                for (i=0;i<n;i++) {
                    nrtk_blsat_t *o;
                    satst_t *st=bl->ss+cs[i].sat-1;
                    if (cs[i].sysidx!=sys||oidx[i]<0||i==r) continue;
                    o=out->sat+oidx[i];
                    if (!o->fixed) continue;
                    for (k=0;k<3;k++) {
                        if (k==2&&(!h3||o->dd_res[2]==0.0)) {o->dd_res[2]=0.0; continue;}
                        o->dd_res[k]+=CLIGHT/cs[i].m->freq[k]*dnr[k];
                    }
                    nrtk_split_iono(cs[i].m->freq[0],cs[i].m->freq[1],o->dd_res[0],o->dd_res[1],
                                    &o->dd_iono,&o->dd_nondisp);
                    for (k=0;k<3;k++) st->res_prev[k]=o->dd_res[k];
                    if (!h3) st->res_has3=0;
                }
                continue;
            }
            nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: reference jump (%d of %d partners)",
                     bl->ista_m,bl->ista_a,cs[r].sat,nj,nt);
            reset_sat(bl->ss+cs[r].sat-1); reset_amb(bl,cs[r].sat);
            bl->nrefix++; bl->nslip++;
            for (i=0;i<n;i++) {
                if (cs[i].sysidx!=sys||oidx[i]<0) continue;
                out->sat[oidx[i]].fixed=out->sat[oidx[i]].isref=out->sat[oidx[i]].newfix=0;
            }
            continue;
        }
        for (i=0;i<n;i++) {
            nrtk_blsat_t *o;
            satst_t *st=bl->ss+cs[i].sat-1;
            if (cs[i].sysidx!=sys||!jump[i]) continue;
            o=out->sat+oidx[i];
            if (jump_integers(cs+i,dj[i],has3[i],dn)) {
                nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d slip repaired (%.0f,%.0f,%.0f)",
                         bl->ista_m,bl->ista_a,cs[i].sat,dn[0],dn[1],dn[2]);
                repair_slip(bl,cs[i].sat,dn,has3[i]);
                bl->nslip++;
                o->slip=1;                      /* the VRS lock time of this satellite restarts */
                for (k=0;k<3;k++) {
                    if (k==2&&(!has3[i]||o->dd_res[2]==0.0)) {o->dd_res[2]=0.0; continue;}
                    o->dd_res[k]-=CLIGHT/cs[i].m->freq[k]*dn[k];
                }
                nrtk_split_iono(cs[i].m->freq[0],cs[i].m->freq[1],o->dd_res[0],o->dd_res[1],
                                &o->dd_iono,&o->dd_nondisp);
                for (k=0;k<3;k++) st->res_prev[k]=o->dd_res[k];
                if (!has3[i]) st->res_has3=0;
                o->ewlfixed=(uint8_t)st->ewl;
                continue;
            }
            nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: jump",bl->ista_m,bl->ista_a,cs[i].sat);
            reset_sat(st); reset_amb(bl,cs[i].sat);
            bl->nrefix++; bl->nslip++;
            o->fixed=o->newfix=0;
        }
    }
    for (i=nfix=0;i<out->ns;i++) nfix+=out->sat[i].fixed;
    out->nsat=n;
    out->nfixed=nfix;
    out->nwl=nwl;
    out->nslip=nslip;
    {
        double pos[3],zhd,zwd;
        ecef2pos(bl->rm,pos); nrtk_trop_zenith(pos,0.7,&zhd,&zwd); out->ztd[0]=zhd+zwd;
        ecef2pos(bl->ra,pos); nrtk_trop_zenith(pos,0.7,&zhd,&zwd); out->ztd[1]=zhd+zwd;
        if (bl->engine==NRTK_ENGINE_RTKLIB&&bl->rtk.opt.tropopt>=TROPOPT_EST) {
            model_zwd(bl->ra,&zwd);
            out->ztd[1]+=bl->rtk.x[IT_(0,&bl->rtk.opt)]-zwd;
        }
        else if (bl->engine==NRTK_ENGINE_NATIVE) out->ztd[1]+=bl->nat.x[0];
    }
    out->stat=nfix>=bl->opt.min_fix+1?2:1;
    update_stat(bl,em->time,out->stat==2);
    return 1;
}
extern void nrtk_bl_unfix(nrtk_baseline_t *bl, int sat)
{
    satst_t *st;
    if (sat<1||sat>MAXSAT) return;
    st=bl->ss+sat-1;
    if (st->nl) bl->nrefix++;
    nrtk_log(NRTK_LOG_DEBUG,"netpos","bl %d-%d sat=%d unfix: network closure",bl->ista_m,
             bl->ista_a,sat);
    st->nl=0; st->nbad=0; st->nrej=0;
    st->wl=0; reset_mw(st,0);
    reset_amb(bl,sat);
}
extern void nrtk_bl_stat(const nrtk_baseline_t *bl, nrtk_blstat_t *st)
{
    int i,nf=0;
    memset(st,0,sizeof(*st));
    st->nepoch=bl->nep;
    st->nfixep=bl->nfixep;
    st->ttff=bl->everfixed?bl->ttff:-1.0;
    st->nrefix=bl->nrefix;
    st->nslip=bl->nslip;
    for (i=0;i<bl->hn;i++) nf+=bl->hist[i];
    st->fixrate_60s=bl->hn>0?(double)nf/bl->hn:0.0;
}
