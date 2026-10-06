/*------------------------------------------------------------------------------
 * atmos.cpp : atmospheric error interpolation (LIM / LSM) and ionosphere
 *             activity monitor (spec ch.6)
 *----------------------------------------------------------------------------*/
#include <cmath>
#include <cstring>
#include <deque>
#include <map>
#include <mutex>
#include <vector>
#include "nrtk/atmos.h"
#include "nrtk/tropo.h"

extern "C" void nrtk_atmos_opt_default(nrtk_atmos_opt_t *o)
{
    memset(o,0,sizeof(*o));
    o->method=NRTK_INTERP_AUTO;
    o->min_bl=2;
    o->trop_model=1;
    o->height_scale=0;
    o->iono_rate_thres=0.1;
    o->max_extrap_km=70.0;
    o->min_sats=6;
}

extern "C" int nrtk_plane_fit(const double *n, const double *e, const double *y,
                              int cnt, double *a, double *b, double *rms)
{
    double nn=0,ne=0,ee=0,ny=0,ey=0,det,ss=0,r;
    int i;

    if (cnt<2) return 0;
    for (i=0;i<cnt;i++) {
        nn+=n[i]*n[i]; ne+=n[i]*e[i]; ee+=e[i]*e[i];
        ny+=n[i]*y[i]; ey+=e[i]*y[i];
    }
    det=nn*ee-ne*ne;
    /* reject near-collinear geometry (relative to the spread of points) */
    if (det<=1E-6*(nn+ee)*(nn+ee)) return 0;
    *a=( ee*ny-ne*ey)/det;
    *b=(-ne*ny+nn*ey)/det;
    for (i=0;i<cnt;i++) {r=y[i]-*a*n[i]-*b*e[i]; ss+=r*r;}
    *rms=cnt>2?sqrt(ss/(cnt-2)):0.0;
    return 1;
}

namespace {

/* point in triangle (plane, km) */
bool in_tri(double px, double py, double x0, double y0, double x1, double y1,
            double x2, double y2)
{
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
    if (fabs(den)<1E-12) return false;
    double l0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
    double l1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
    double l2=1.0-l0-l1;
    return l0>=-1E-6&&l1>=-1E-6&&l2>=-1E-6;
}
/* V inside the hull of master (origin) and the aux points used */
bool in_hull(double vn, double ve, const std::vector<double> &n, const std::vector<double> &e)
{
    std::vector<double> x(1,0.0), y(1,0.0);
    x.insert(x.end(),e.begin(),e.end()); y.insert(y.end(),n.begin(),n.end());
    size_t m=x.size();
    for (size_t i=0;i<m;i++) for (size_t j=i+1;j<m;j++) for (size_t k=j+1;k<m;k++) {
        if (in_tri(ve,vn,x[i],y[i],x[j],y[j],x[k],y[k])) return true;
    }
    return false;
}
const nrtk_blsat_t *find_sat(const nrtk_blepoch_t *bl, int sat)
{
    for (int i=0;i<bl->ns;i++) if (bl->sat[i].sat==sat) return bl->sat+i;
    return nullptr;
}

} // namespace

#define TROP_SCALE_H   2000.0   /* scale height of the wet delay (m) */
#define ZWD_ANOM_PRIOR 0.30     /* prior std of the master's zenith wet anomaly (m): weak, it only
                                   keeps c at zero when the cell has no height spread */
#define ND_SIGMA       0.015    /* typical non-dispersive DD residual (m), scales the prior */
#define ZWD_ANOM_RW    0.02     /* random walk of the sea-level anomaly (m/sqrt(h)) */
#define ZWD_ANOM_MAXDT 7200.0   /* state older than this is started over (s) */

namespace {
/* one satellite of the cell: DD values master->aux_i and their partials */
struct satrow {
    int k;                      /* index in the master epoch */
    int sys;
    int isref,refcnt,refslip;
    int usable;                 /* enough baselines */
    double gm;                  /* wet mapping function difference sat - ref at the master */
    std::vector<double> n,e,I,R,G;
    uint32_t mask;
    int newfix,slipped;
};
}

extern "C" int nrtk_atmos_interp(const nrtk_atmos_opt_t *opt,
                                 const nrtk_staepoch_t *mep,
                                 const nrtk_cellbl_t *bls, int nbl,
                                 const double *vpos, nrtk_zwd_state_t *st,
                                 nrtk_vcorr_t *out)
{
    const double ep2000[]={2000,1,1,12,0,0};
    double vn,ve,posv[3],posm[3],ss_nd=0.0,mjd,fv=1.0,mh,mw,c=0.0;
    int nnd=0;
    std::vector<double> bn(nbl),be(nbl),bf(nbl,1.0),mwm(mep->ns,0.0);
    std::vector<satrow> rows;

    memset(out,0,sizeof(*out));
    out->time=mep->time;
    out->master=mep->sta;
    if (nbl<1) return 0;

    nrtk_ne(mep->rr,vpos,&vn,&ve); vn/=1E3; ve/=1E3;
    ecef2pos(vpos,posv); ecef2pos(mep->rr,posm);
    /* height_scale (experimental, off by default): the wet delay anomaly (true - model)
       of a station is assumed to decay with its height like the wet delay itself,
       anomaly(h) = anomaly(h_master) * f(h), f = exp(-(h-h_m)/H). A plane through the
       lowland stations alone puts the lowland anomaly on a mountain VRS (Taiwan BAL2,
       1699 m above its cell: 16 cm height error of the rover, reduced to 3 cm with this
       model), but a single mountain station in cloud sets a wrong profile for the whole
       cell (northern ring: ZUZH/ANBU residuals 1.8 -> 3.5 cm), see VALIDATION 9.4. */
    for (int i=0;i<nbl;i++) {
        double pa[3];
        nrtk_ne(mep->rr,bls[i].apos,&bn[i],&be[i]);
        bn[i]/=1E3; be[i]/=1E3;
        ecef2pos(bls[i].apos,pa);
        bf[i]=exp(-(pa[2]-posm[2])/TROP_SCALE_H);
    }
    fv=exp(-(posv[2]-posm[2])/TROP_SCALE_H);
    mjd=51544.5+timediff(mep->time,epoch2time(ep2000))/86400.0;
    if (opt->trop_model) {
        for (int k=0;k<mep->ns;k++) {
            nrtk_gmf(mjd,posm[0],posm[1],posm[2],PI/2.0-mep->sat[k].azel[1],&mh,&mw);
            mwm[k]=mw;
        }
    }

    /* pass 1: reference satellite per system, DD values of the others */
    for (int sys=0;sys<NRTK_NSYS;sys++) {
        /* common reference: satellite fixed in most baselines, then elevation */
        int ref=0,refcnt=0,kref=-1; double refel=-1.0;
        for (int k=0;k<mep->ns;k++) {
            const nrtk_satdata_t *s=mep->sat+k;
            if (s->sysidx!=sys) continue;
            int cnt=0;
            for (int i=0;i<nbl;i++) {
                if (!bls[i].bl||!bls[i].bl->valid) continue;
                const nrtk_blsat_t *b=find_sat(bls[i].bl,s->sat);
                if (b&&b->fixed) cnt++;
            }
            if (cnt>refcnt||(cnt==refcnt&&cnt>0&&s->azel[1]>refel)) {
                ref=s->sat; refcnt=cnt; refel=s->azel[1]; kref=k;
            }
        }
        out->refsat[sys]=ref;
        if (!ref||refcnt<opt->min_bl) continue;

        for (int k=0;k<mep->ns;k++) {
            const nrtk_satdata_t *s=mep->sat+k;
            if (s->sysidx!=sys) continue;
            satrow row;
            row.k=k; row.sys=sys; row.isref=(s->sat==ref); row.refcnt=refcnt; row.refslip=0;
            row.usable=0; row.gm=0.0; row.mask=0; row.newfix=row.slipped=0;
            if (row.isref) {
                for (int i=0;i<nbl;i++) {
                    if (!bls[i].bl||!bls[i].bl->valid) continue;
                    const nrtk_blsat_t *b=find_sat(bls[i].bl,ref);
                    if (b&&b->fixed&&b->slip) row.refslip=1;
                }
                rows.push_back(row);
                continue;
            }
            row.gm=opt->trop_model&&kref>=0?mwm[k]-mwm[kref]:0.0;
            /* DD values master->aux_i w.r.t. the common reference */
            for (int i=0;i<nbl;i++) {
                const nrtk_blepoch_t *bl=bls[i].bl;
                if (!bl||!bl->valid) continue;
                const nrtk_blsat_t *bs=find_sat(bl,s->sat),*br=find_sat(bl,ref);
                if (!bs||!br||!bs->fixed||!br->fixed) continue;
                double sg=bls[i].flip?-1.0:1.0;
                double ddi=sg*(bs->dd_iono-br->dd_iono);
                double ddd=sg*(bs->dd_nondisp-br->dd_nondisp);
                double ddt=sg*(bs->dd_trop_model-br->dd_trop_model);
                row.n.push_back(bn[i]); row.e.push_back(be[i]);
                row.I.push_back(ddi);
                row.R.push_back(opt->trop_model?ddd-ddt:ddd);
                row.G.push_back((bf[i]-1.0)*row.gm);
                row.mask|=1u<<i;
                row.newfix|=bs->newfix|br->newfix;
                row.slipped|=bs->slip;
            }
            int cnt=(int)row.n.size();
            row.usable=cnt>=opt->min_bl;
            if (row.usable&&opt->method==NRTK_INTERP_LIM&&cnt>2) {   /* triangle partners only */
                row.n.resize(2); row.e.resize(2); row.I.resize(2); row.R.resize(2); row.G.resize(2);
                uint32_t m2=0; int cc=0;
                for (int i=0;i<nbl&&cc<2;i++) if (row.mask&(1u<<i)) {m2|=1u<<i; cc++;}
                row.mask=m2;
            }
            rows.push_back(row);
        }
    }

    /* master zenith wet anomaly c from the height spread of the cell: for every
       satellite R_i = a*n_i + b*e_i + c*(f_i-1)*(mw_s-mw_r); the plane coefficients
       are eliminated (Schur complement), c is common to all satellites. Without
       height differences the prior keeps c at zero (= the plain plane). */
    if (opt->trop_model&&opt->height_scale) {
        double Ncc=(ND_SIGMA/ZWD_ANOM_PRIOR)*(ND_SIGMA/ZWD_ANOM_PRIOR),bc=0.0,Ninfo=0.0;
        for (const satrow &r:rows) {
            if (r.isref||!r.usable) continue;
            int cnt=(int)r.n.size();
            double nn=0,ne=0,ee=0,ng=0,eg=0,nR=0,eR=0,gg=0,gR=0;
            for (int i=0;i<cnt;i++) {
                nn+=r.n[i]*r.n[i]; ne+=r.n[i]*r.e[i]; ee+=r.e[i]*r.e[i];
                ng+=r.n[i]*r.G[i]; eg+=r.e[i]*r.G[i];
                nR+=r.n[i]*r.R[i]; eR+=r.e[i]*r.R[i];
                gg+=r.G[i]*r.G[i]; gR+=r.G[i]*r.R[i];
            }
            double det=nn*ee-ne*ne;
            if (det<=1E-6*(nn+ee)*(nn+ee)) continue;
            double i11=ee/det,i12=-ne/det,i22=nn/det;
            double un=i11*ng+i12*eg,ue=i12*ng+i22*eg;      /* N^-1 u */
            Ncc+=gg-(ng*un+eg*ue);
            bc+=gR-(nR*un+eR*ue);
        }
        Ninfo=Ncc-(ND_SIGMA/ZWD_ANOM_PRIOR)*(ND_SIGMA/ZWD_ANOM_PRIOR);
        c=bc/Ncc;
        if (fabs(c)>1.0) {c=0.0; Ninfo=0.0;}                /* implausible */
        if (st) {
            /* scalar Kalman filter on the sea-level anomaly c0 = c*exp(h_m/H): the
               per-epoch estimate (variance ND_SIGMA^2/Ninfo) is noisy when the cell
               has little height spread, the anomaly itself changes slowly */
            double sc=exp(posm[2]/TROP_SCALE_H),dt;
            if (st->valid&&(dt=timediff(mep->time,st->time))>=0.0&&dt<=ZWD_ANOM_MAXDT) {
                st->var+=ZWD_ANOM_RW*ZWD_ANOM_RW/3600.0*dt;
            }
            else {st->c0=0.0; st->var=ZWD_ANOM_PRIOR*ZWD_ANOM_PRIOR; st->valid=1;}
            if (Ninfo>1E-9) {
                double R=ND_SIGMA*ND_SIGMA/Ninfo*sc*sc,K=st->var/(st->var+R);
                st->c0+=K*(c*sc-st->c0);
                st->var*=1.0-K;
            }
            st->time=mep->time;
            c=st->c0/sc;
        }
    }
    out->zwd_anom=c;

    /* pass 2: planes per satellite, evaluated at the VRS */
    for (satrow &r:rows) {
        if (out->ns>=NRTK_BL_MAXSAT) break;
        const nrtk_satdata_t *s=mep->sat+r.k;
        nrtk_vsat_t *o=out->sat+out->ns++;
        memset(o,0,sizeof(*o));
        o->sat=s->sat; o->sysidx=(int8_t)r.sys;
        if (r.isref) {
            o->ok=1; o->isref=1; o->nbl=(uint8_t)r.refcnt; o->slip=(uint8_t)r.refslip;
            out->nok++;
            continue;
        }
        if (!r.usable) continue;
        int cnt=(int)r.n.size();
        for (int i=0;i<cnt;i++) r.R[i]-=c*r.G[i];
        double ai,bi,ri,an,bnn,rn;
        if (!nrtk_plane_fit(r.n.data(),r.e.data(),r.I.data(),cnt,&ai,&bi,&ri)) continue;
        if (!nrtk_plane_fit(r.n.data(),r.e.data(),r.R.data(),cnt,&an,&bnn,&rn)) continue;

        /* model troposphere master->VRS */
        double ddt_v=0.0;
        if (opt->trop_model) {
            const nrtk_satdata_t *sr=nrtk_ep_sat(mep,out->refsat[r.sys]);
            double azv[2],azr[2],ev[3],er[3],dv=0,dr=0;
            for (int j=0;j<3;j++) {ev[j]=s->rs[j]-vpos[j]; er[j]=sr->rs[j]-vpos[j];}
            dv=norm(ev,3); dr=norm(er,3);
            for (int j=0;j<3;j++) {ev[j]/=dv; er[j]/=dr;}
            satazel(posv,ev,azv); satazel(posv,er,azr);
            ddt_v=(nrtk_trop_slant(mep->time,posv,azv)-nrtk_trop_slant(mep->time,posm,s->azel))-
                  (nrtk_trop_slant(mep->time,posv,azr)-nrtk_trop_slant(mep->time,posm,sr->azel));
        }
        o->ok=1;
        o->nbl=(uint8_t)cnt;
        o->blmask=r.mask;
        o->newfix=(uint8_t)r.newfix;
        o->slip=(uint8_t)r.slipped;
        o->dd_iono=ai*vn+bi*ve;
        o->dd_trop_model=ddt_v;
        o->dd_nondisp=an*vn+bnn*ve+c*(fv-1.0)*r.gm+ddt_v;
        o->rms_iono=ri; o->rms_nd=rn;
        o->coef_iono[0]=ai; o->coef_iono[1]=bi;
        o->coef_nd[0]=an; o->coef_nd[1]=bnn;
        o->extrap=(uint8_t)!in_hull(vn,ve,r.n,r.e);
        for (double rr:r.R) {ss_nd+=rr*rr; nnd++;}
        out->nok++;
    }
    out->trop_resid=nnd>0?sqrt(ss_nd/nnd):0.0;
    /* per-satellite extrapolation is a quality indicator only; the cell level
       extrapolation flag is set by the caller (cell outside the network) */
    if (out->nok<opt->min_sats) out->degraded|=NRTK_DEG_FEWSAT;
    return out->nok;
}

/* ionosphere activity monitor --------------------------------------------------*/
struct nrtk_ionomon {
    double window;
    struct Sample { gtime_t t; double I; int ref; };
    struct BlState {
        std::map<int,std::deque<Sample>> sat;
        double rate=0.0;
        gtime_t t{};
    };
    std::map<std::pair<int,int>,BlState> bl;
    mutable std::mutex mtx;
};

extern "C" nrtk_ionomon_t *nrtk_ionomon_create(double window_s)
{
    nrtk_ionomon_t *m=new nrtk_ionomon_t;
    m->window=window_s>0.0?window_s:60.0;
    return m;
}
extern "C" void nrtk_ionomon_free(nrtk_ionomon_t *m) { delete m; }

extern "C" void nrtk_ionomon_add(nrtk_ionomon_t *m, const nrtk_blepoch_t *b)
{
    if (!b->valid) return;
    std::lock_guard<std::mutex> lk(m->mtx);
    auto key=std::make_pair(std::min(b->sta_m,b->sta_a),std::max(b->sta_m,b->sta_a));
    auto &st=m->bl[key];
    double rmax=0.0;
    for (int i=0;i<b->ns;i++) {
        const nrtk_blsat_t *s=b->sat+i;
        int ref=b->refsat[s->sysidx];
        auto &q=st.sat[s->sat];
        if (!s->fixed||s->isref||!ref) {q.clear(); continue;}
        /* restart after a reference change or a gap longer than the window (the
           rate is taken between the oldest and the newest sample, so the
           sampling interval itself does not matter: 1 s to 30 s data) */
        if (!q.empty()&&(q.back().ref!=ref||timediff(b->time,q.back().t)>m->window)) q.clear();
        q.push_back({b->time,s->dd_iono,ref});
        while (q.size()>2&&timediff(b->time,q.front().t)>m->window) q.pop_front();
        double dt=timediff(b->time,q.front().t);
        if (dt>=0.8*m->window&&s->azel_m[1]>=20.0*D2R) {
            double r=fabs(s->dd_iono-q.front().I)/dt*60.0;
            if (r>rmax) rmax=r;
        }
    }
    st.rate=rmax;
    st.t=b->time;
}
extern "C" double nrtk_ionomon_rate(const nrtk_ionomon_t *m, int sta_m, int sta_a)
{
    std::lock_guard<std::mutex> lk(m->mtx);
    auto key=std::make_pair(std::min(sta_m,sta_a),std::max(sta_m,sta_a));
    auto it=m->bl.find(key);
    return it==m->bl.end()?0.0:it->second.rate;
}
