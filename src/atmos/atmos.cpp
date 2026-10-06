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

extern "C" int nrtk_atmos_interp(const nrtk_atmos_opt_t *opt,
                                 const nrtk_staepoch_t *mep,
                                 const nrtk_cellbl_t *bls, int nbl,
                                 const double *vpos, nrtk_vcorr_t *out)
{
    double vn,ve,posv[3],posm[3],ss_nd=0.0;
    int nnd=0;
    std::vector<double> bn(nbl),be(nbl);

    memset(out,0,sizeof(*out));
    out->time=mep->time;
    out->master=mep->sta;
    if (nbl<1) return 0;

    nrtk_ne(mep->rr,vpos,&vn,&ve); vn/=1E3; ve/=1E3;
    for (int i=0;i<nbl;i++) {
        nrtk_ne(mep->rr,bls[i].apos,&bn[i],&be[i]);
        bn[i]/=1E3; be[i]/=1E3;
    }
    ecef2pos(vpos,posv); ecef2pos(mep->rr,posm);

    for (int sys=0;sys<NRTK_NSYS;sys++) {
        /* common reference: satellite fixed in most baselines, then elevation */
        int ref=0,refcnt=0; double refel=-1.0;
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
                ref=s->sat; refcnt=cnt; refel=s->azel[1];
            }
        }
        out->refsat[sys]=ref;
        if (!ref||refcnt<opt->min_bl) continue;

        for (int k=0;k<mep->ns&&out->ns<NRTK_BL_MAXSAT;k++) {
            const nrtk_satdata_t *s=mep->sat+k;
            if (s->sysidx!=sys) continue;
            nrtk_vsat_t *o=out->sat+out->ns++;
            memset(o,0,sizeof(*o));
            o->sat=s->sat; o->sysidx=(int8_t)sys;
            if (s->sat==ref) {
                o->ok=1; o->isref=1; o->nbl=(uint8_t)refcnt; out->nok++;
                for (int i=0;i<nbl;i++) {
                    if (!bls[i].bl||!bls[i].bl->valid) continue;
                    const nrtk_blsat_t *b=find_sat(bls[i].bl,ref);
                    if (b&&b->fixed&&b->slip) o->slip=1;
                }
                continue;
            }

            /* DD values master->aux_i w.r.t. the common reference */
            std::vector<double> n,e,I,R;
            uint32_t mask=0;
            int newfix=0,slipped=0;
            for (int i=0;i<nbl;i++) {
                const nrtk_blepoch_t *bl=bls[i].bl;
                if (!bl||!bl->valid) continue;
                const nrtk_blsat_t *bs=find_sat(bl,s->sat),*br=find_sat(bl,ref);
                if (!bs||!br||!bs->fixed||!br->fixed) continue;
                double sg=bls[i].flip?-1.0:1.0;
                double ddi=sg*(bs->dd_iono-br->dd_iono);
                double ddd=sg*(bs->dd_nondisp-br->dd_nondisp);
                double ddt=sg*(bs->dd_trop_model-br->dd_trop_model);
                n.push_back(bn[i]); e.push_back(be[i]);
                I.push_back(ddi);
                R.push_back(opt->trop_model?ddd-ddt:ddd);
                mask|=1u<<i;
                newfix|=bs->newfix|br->newfix;
                slipped|=bs->slip;
            }
            int cnt=(int)n.size();
            if (cnt<opt->min_bl) continue;
            if (opt->method==NRTK_INTERP_LIM&&cnt>2) {   /* triangle partners only */
                n.resize(2); e.resize(2); I.resize(2); R.resize(2); cnt=2;
                uint32_t m2=0; int c=0;
                for (int i=0;i<nbl&&c<2;i++) if (mask&(1u<<i)) {m2|=1u<<i; c++;}
                mask=m2;
            }
            double ai,bi,ri,an,bnn,rn;
            if (!nrtk_plane_fit(n.data(),e.data(),I.data(),cnt,&ai,&bi,&ri)) continue;
            if (!nrtk_plane_fit(n.data(),e.data(),R.data(),cnt,&an,&bnn,&rn)) continue;

            /* model troposphere master->VRS */
            double ddt_v=0.0;
            if (opt->trop_model) {
                const nrtk_satdata_t *sr=nrtk_ep_sat(mep,ref);
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
            o->blmask=mask;
            o->newfix=(uint8_t)newfix;
            o->slip=(uint8_t)slipped;
            o->dd_iono=ai*vn+bi*ve;
            o->dd_trop_model=ddt_v;
            o->dd_nondisp=an*vn+bnn*ve+ddt_v;
            o->rms_iono=ri; o->rms_nd=rn;
            o->coef_iono[0]=ai; o->coef_iono[1]=bi;
            o->coef_nd[0]=an; o->coef_nd[1]=bnn;
            o->extrap=(uint8_t)!in_hull(vn,ve,n,e);
            for (double r:R) {ss_nd+=r*r; nnd++;}
            out->nok++;
        }
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
        /* restart after a reference change or a gap (half the window: the rate
           is taken between the oldest and the newest sample, so the sampling
           interval itself does not matter) */
        if (!q.empty()&&(q.back().ref!=ref||timediff(b->time,q.back().t)>0.5*m->window)) q.clear();
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
