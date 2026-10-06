/*------------------------------------------------------------------------------
 * network.cpp : processing pipeline shared by offline and real-time modes
 *----------------------------------------------------------------------------*/
#include <algorithm>
#include <cmath>
#include <tuple>
#include <chrono>
#include <cstring>
#include "nrtk/pipeline.hpp"

namespace nrtk {

/* thread pool -------------------------------------------------------------------*/
ThreadPool::ThreadPool(int n)
{
    if (n<1) n=1;
    for (int i=0;i<n;i++) th_.emplace_back([this]{worker();});
}
ThreadPool::~ThreadPool()
{
    {std::lock_guard<std::mutex> lk(mtx_); stop_=true;}
    cv_.notify_all();
    for (auto &t:th_) t.join();
}
void ThreadPool::worker()
{
    long seen=0;
    for (;;) {
        std::function<void(int)> job;
        int i;
        {
            std::unique_lock<std::mutex> lk(mtx_);
            cv_.wait(lk,[&]{return stop_||(gen_!=seen&&next_<total_);});
            if (stop_) return;
            i=next_++;
            if (next_>=total_) seen=gen_;
            job=job_;
        }
        job(i);
        {
            std::lock_guard<std::mutex> lk(mtx_);
            if (--pending_==0) done_cv_.notify_all();
        }
    }
}
void ThreadPool::parallel_for(int n, const std::function<void(int)> &f)
{
    if (n<=0) return;
    if (th_.size()<=1||n==1) {for (int i=0;i<n;i++) f(i); return;}
    std::unique_lock<std::mutex> lk(mtx_);
    job_=f; next_=0; total_=n; pending_=n; gen_++;
    cv_.notify_all();
    done_cv_.wait(lk,[&]{return pending_==0;});
}

std::string time_str(gtime_t t, int ndec)
{
    char s[64];
    time2str(t,s,ndec);
    return s;
}

const nrtk_blepoch_t *NetEpoch::find(int a, int b, int *flip) const
{
    auto it=bl.find({a,b});
    if (it!=bl.end()) {if (flip) *flip=0; return &it->second;}
    it=bl.find({b,a});
    if (it!=bl.end()) {if (flip) *flip=1; return &it->second;}
    return nullptr;
}

/* network ------------------------------------------------------------------------*/
Network::Network(const nrtk_config_t &cfg, const nrtk_atx_t *atx) : cfg_(cfg), atx_(atx)
{
    ant_.resize(cfg_.nsta,nullptr);
    st_.resize(cfg_.nsta);
    last_data_.resize(cfg_.nsta,gtime_t{});
    for (int i=0;i<cfg_.nsta;i++) {
        st_[i].id=cfg_.sta[i].id;
        if (atx_&&cfg_.sta[i].antenna[0]) {
            ant_[i]=nrtk_atx_find(atx_,cfg_.sta[i].antenna);
            if (!ant_[i]) {
                nrtk_log(NRTK_LOG_WARN,"network","antenna '%s' of %s not in ANTEX: PCO/PCV not applied",
                         cfg_.sta[i].antenna,cfg_.sta[i].id);
            }
        }
    }
    cm_=nrtk_cellmgr_create(&cfg_.cell,cfg_.sta,cfg_.nsta);
    int (*pair)[2]=new int[NRTK_MAXSTA*NRTK_MAXSTA][2];
    int np=nrtk_cellmgr_baselines(cm_,pair,NRTK_MAXSTA*NRTK_MAXSTA);
    for (int k=0;k<np;k++) {
        int i=pair[k][0],j=pair[k][1];
        nrtk_baseline_t *b=nrtk_bl_create(&cfg_.netpos,cfg_.sta+i,i,cfg_.sta+j,j);
        if (b) bls_.push_back({i,j,b});
    }
    delete[] pair;
    mon_=nrtk_ionomon_create(60.0);
    int nth=cfg_.threads>0?cfg_.threads:(int)std::thread::hardware_concurrency();
    if (nth<1) nth=1;
    if (nth>16) nth=16;
    pool_.reset(new ThreadPool(nth));
    nrtk_log(NRTK_LOG_INFO,"network","%d stations, %zu baselines, %d threads",
             cfg_.nsta,bls_.size(),nth);
}
Network::~Network()
{
    pool_.reset();
    for (auto &b:bls_) nrtk_bl_free(b.bl);
    nrtk_cellmgr_free(cm_);
    nrtk_ionomon_free(mon_);
}

std::shared_ptr<NetEpoch> Network::process(gtime_t t,
                                           const std::vector<std::vector<obsd_t>> &raw,
                                           std::shared_ptr<const nav_t> navp,
                                           const std::vector<double> *latency)
{
    const nav_t &nav=*navp;
    auto t0=std::chrono::steady_clock::now();
    auto ep=std::make_shared<NetEpoch>();
    ep->t=t;
    ep->nav=navp;
    ep->sta.resize(cfg_.nsta);
    ep->latency.assign(cfg_.nsta,0.0);
    if (latency) ep->latency=*latency;
    if (!tstart_.time) tstart_=t;

    /* prepare station epochs */
    pool_->parallel_for(cfg_.nsta,[&](int i){
        if (i>=(int)raw.size()||raw[i].empty()||!cfg_.sta[i].enabled) return;
        auto se=std::make_shared<nrtk_staepoch_t>();
        if (nrtk_prep_epoch(&cfg_.prep,cfg_.sta+i,i,ant_[i],raw[i].data(),
                            (int)raw[i].size(),&nav,se.get())>=4) {
            ep->sta[i]=se;
        }
    });
    /* baselines */
    std::vector<nrtk_blepoch_t> res(bls_.size());
    std::vector<int> ok(bls_.size(),0);
    pool_->parallel_for((int)bls_.size(),[&](int k){
        const BL &b=bls_[k];
        const nrtk_staepoch_t *em=ep->sta[b.m].get(),*ea=ep->sta[b.a].get();
        if (!em||!ea) {
            /* keep statistics running (epoch without data) */
            memset(&res[k],0,sizeof(res[k]));
            return;
        }
        ok[k]=nrtk_bl_process(b.bl,em,ea,&nav,&res[k]);
    });
    for (size_t k=0;k<bls_.size();k++) {
        if (!ok[k]) continue;
        ep->bl[{bls_[k].m,bls_[k].a}]=res[k];
    }
    if (cfg_.netpos.closure) closure_check(*ep);
    for (auto &kv:ep->bl) nrtk_ionomon_add(mon_,&kv.second);
    update_health(*ep);
    ep->proc_ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t0).count();
    return ep;
}

/* network closure of the fixed DD residuals --------------------------------------
 * Around a triangle i->j->k->i every observation enters twice with opposite
 * sign, so the closure of correctly fixed DD phase residuals is exactly zero;
 * a wrong integer leaves a multiple of the wavelength. With each baseline's
 * own reference satellite the closure is the same constant for all
 * consistent satellites of a system; outliers from the median are wrong.
 * A (baseline, satellite) fix that fails in every triangle it takes part in
 * (and passes in none) is dropped and resolved again. */
static const nrtk_blsat_t *blsat(const nrtk_blepoch_t *b, int sat)
{
    for (int i=0;i<b->ns;i++) if (b->sat[i].sat==sat) return b->sat+i;
    return nullptr;
}
void Network::closure_check(NetEpoch &ep)
{
    const double tol=0.01;
    int n=cfg_.nsta;
    struct Vote {int bad=0,good=0;};
    std::map<std::tuple<int,int,int>,Vote> votes;   /* (stored m, stored a, sat) */

    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) for (int k=j+1;k<n;k++) {
        int fl[3];
        const nrtk_blepoch_t *b[3]={ep.find(i,j,fl),ep.find(j,k,fl+1),ep.find(k,i,fl+2)};
        if (!b[0]||!b[1]||!b[2]) continue;
        for (int sys=0;sys<NRTK_NSYS;sys++) {
            std::vector<int> sats;
            std::vector<double> c0,c1;
            for (int q=0;q<b[0]->ns;q++) {
                const nrtk_blsat_t *s0=b[0]->sat+q;
                if (s0->sysidx!=sys||!s0->fixed) continue;
                const nrtk_blsat_t *s1=blsat(b[1],s0->sat),*s2=blsat(b[2],s0->sat);
                if (!s1||!s2||!s1->fixed||!s2->fixed) continue;
                const nrtk_blsat_t *ss[3]={s0,s1,s2};
                double x0=0,x1=0;
                for (int e=0;e<3;e++) {
                    double sg=fl[e]?-1.0:1.0;
                    x0+=sg*ss[e]->dd_res[0]; x1+=sg*ss[e]->dd_res[1];
                }
                sats.push_back(s0->sat); c0.push_back(x0); c1.push_back(x1);
            }
            if (sats.size()<3) continue;
            std::vector<double> m0=c0,m1=c1;
            std::nth_element(m0.begin(),m0.begin()+m0.size()/2,m0.end());
            std::nth_element(m1.begin(),m1.begin()+m1.size()/2,m1.end());
            double med0=m0[m0.size()/2],med1=m1[m1.size()/2];
            for (size_t q=0;q<sats.size();q++) {
                bool bad=fabs(c0[q]-med0)>tol||fabs(c1[q]-med1)>tol;
                if (bad&&nrtk_log_enabled(NRTK_LOG_DEBUG)) {
                    char id[8]; satno2id(sats[q],id);
                    nrtk_log(NRTK_LOG_DEBUG,"closure","%s %s %s-%s-%s L1 %.4f L2 %.4f m",
                             time_str(ep.t).c_str(),id,cfg_.sta[i].id,cfg_.sta[j].id,cfg_.sta[k].id,
                             c0[q]-med0,c1[q]-med1);
                }
                for (int e=0;e<3;e++) {
                    Vote &v=votes[{b[e]->sta_m,b[e]->sta_a,sats[q]}];
                    if (bad) v.bad++; else v.good++;
                }
            }
        }
    }
    for (auto &kv:votes) {
        if (kv.second.bad==0||kv.second.good>0) continue;
        int m=std::get<0>(kv.first),a=std::get<1>(kv.first),sat=std::get<2>(kv.first);
        auto it=ep.bl.find({m,a});
        if (it==ep.bl.end()) continue;
        for (int q=0;q<it->second.ns;q++) {
            nrtk_blsat_t &s=it->second.sat[q];
            if (s.sat!=sat) continue;
            s.fixed=0;
            it->second.nfixed--;
        }
        for (auto &b:bls_) if (b.m==m&&b.a==a) nrtk_bl_unfix(b.bl,sat);
        nclosure_++;
        char id[8]; satno2id(sat,id);
        nrtk_log(NRTK_LOG_INFO,"network","%s closure: %s-%s %s fix dropped (fails %d triangles)",
                 time_str(ep.t).c_str(),cfg_.sta[m].id,cfg_.sta[a].id,id,kv.second.bad);
    }
}

void Network::update_health(NetEpoch &ep)
{
    std::lock_guard<std::mutex> lk(mtx_);
    ep.healthy.assign(cfg_.nsta,0);
    bool warmup=timediff(ep.t,tstart_)<180.0;

    for (int i=0;i<cfg_.nsta;i++) {
        StationStatus &s=st_[i];
        bool data=ep.sta[i]!=nullptr;
        if (data) last_data_[i]=ep.t;
        s.last_data_age=last_data_[i].time?timediff(ep.t,last_data_[i]):1E9;
        s.data=data;
        s.latency=ep.latency[i];
        /* best 60 s fix rate among the station's baselines */
        double fr=0.0;
        for (auto &b:bls_) {
            if (b.m!=i&&b.a!=i) continue;
            nrtk_blstat_t bs;
            nrtk_bl_stat(b.bl,&bs);
            fr=std::max(fr,bs.fixrate_60s);
        }
        s.fixrate=fr;
        bool offline=s.last_data_age>cfg_.station_offline_s;
        bool late=s.latency>cfg_.health_latency_s;
        if (offline) s.recovering=true;
        bool healthy;
        if (offline||late||!cfg_.sta[i].enabled) healthy=false;
        else if (s.recovering) {
            healthy=fr>=cfg_.health_fixrate;
            if (healthy) s.recovering=false;
        }
        else if (warmup) healthy=true;
        else if (s.healthy) healthy=fr>=cfg_.health_fixrate-0.15;  /* hysteresis */
        else healthy=fr>=cfg_.health_fixrate;
        if (healthy!=s.healthy) {
            nrtk_log(NRTK_LOG_INFO,"health","%s station %s %s (fixrate %.2f, age %.1f s)",
                     time_str(ep.t).c_str(),s.id.c_str(),healthy?"healthy":"unhealthy",fr,
                     s.last_data_age);
        }
        s.healthy=healthy;
        ep.healthy[i]=healthy;
        nrtk_cellmgr_set_health(cm_,i,healthy);
    }
}

std::vector<StationStatus> Network::station_status() const
{
    std::lock_guard<std::mutex> lk(mtx_);
    return st_;
}
std::vector<BaselineStatus> Network::baseline_status() const
{
    std::vector<BaselineStatus> v;
    for (auto &b:bls_) {
        nrtk_blstat_t s;
        nrtk_bl_stat(b.bl,&s);
        BaselineStatus o;
        o.a=b.m; o.b=b.a;
        o.length_km=nrtk_dist(cfg_.sta[b.m].pos,cfg_.sta[b.a].pos)/1E3;
        o.fixrate_60s=s.fixrate_60s; o.ttff=s.ttff;
        o.nepoch=s.nepoch; o.nfixep=s.nfixep; o.nslip=s.nslip; o.nrefix=s.nrefix;
        o.iono_rate=nrtk_ionomon_rate(mon_,b.m,b.a);
        v.push_back(o);
    }
    return v;
}

/* VRS stream ----------------------------------------------------------------------*/
VrsStream::VrsStream(Network &net, const double *pos, int exclude) : net_(net), exclude_(exclude)
{
    memcpy(pos_,pos,sizeof(pos_));
    const nrtk_config_t &cfg=net.config();
    std::vector<nrtk_station_t> sta(cfg.sta,cfg.sta+cfg.nsta);
    if (exclude>=0&&exclude<cfg.nsta) sta[exclude].enabled=0;
    cm_=nrtk_cellmgr_create(&cfg.cell,sta.data(),cfg.nsta);
    nrtk_vrs_opt_t vo=cfg.vrs;
    memcpy(vo.pos,pos,sizeof(vo.pos));
    vrs_=nrtk_vrs_create(&vo);
    health_.assign(cfg.nsta,1);
}
VrsStream::~VrsStream()
{
    nrtk_cellmgr_free(cm_);
    nrtk_vrs_free(vrs_);
}
void VrsStream::reset() { nrtk_vrs_reset(vrs_); }

/* corrected satellites below which the VRS falls back to the nearest station,
 * and needed to resume the VRS */
static const int VRS_MIN_SATS=5, VRS_RESUME_SATS=8;

int VrsStream::generate(const NetEpoch &ep, VrsOutput &out)
{
    const nrtk_config_t &cfg=net_.config();
    out.n=0; out.mode=0; out.degraded=0; out.iono_activity=0.0;
    memset(&out.corr,0,sizeof(out.corr));
    memset(&out.info,0,sizeof(out.info));

    for (int i=0;i<cfg.nsta;i++) {
        int h=i!=exclude_&&i<(int)ep.healthy.size()&&ep.healthy[i];
        if (h!=health_[i]) {nrtk_cellmgr_set_health(cm_,i,h); health_[i]=h;}
    }
    nrtk_cellmgr_select(cm_,pos_,&out.cell);
    nrtk_cell_t &cell=out.cell;

    if (cell.valid) {
        /* master must have data in this epoch, else another triangle vertex */
        int master=cell.master;
        if (!ep.sta[master]) {
            master=-1;
            for (int k=0;k<3;k++) {
                int v=cell.tri[k];
                if (v>=0&&v!=cell.master&&ep.sta[v]) {master=v; break;}
            }
            if (master<0) cell.valid=0;
            else { /* swap roles */
                for (int k=0;k<cell.naux;k++) if (cell.aux[k]==master) cell.aux[k]=cell.master;
                cell.master=master;
                cell.dist_km=nrtk_dist(pos_,cfg.sta[master].pos)/1E3;
            }
        }
    }
    if (cell.valid) {
        std::vector<nrtk_cellbl_t> cb;
        for (int k=0;k<cell.naux;k++) {
            int a=cell.aux[k],flip=0;
            const nrtk_blepoch_t *b=ep.find(cell.master,a,&flip);
            if (!b||!ep.sta[a]) continue;
            nrtk_cellbl_t c;
            c.bl=b;
            /* stored key (m,a): DD values are aux-side minus master-side of
               the stored baseline; flip when the stored master side is aux */
            c.flip=b->sta_m!=cell.master;
            c.aux=a;
            nrtk_sta_pos(cfg.sta+a,ep.t,c.apos);
            cb.push_back(c);
            out.iono_activity=std::max(out.iono_activity,nrtk_ionomon_rate(net_.ionomon(),cell.master,a));
        }
        const nrtk_staepoch_t *mep=ep.sta[cell.master].get();
        nrtk_atmos_interp(&cfg.atmos,mep,cb.data(),(int)cb.size(),pos_,&out.corr);
        out.corr.iono_activity=out.iono_activity;
        if (out.iono_activity>cfg.atmos.iono_rate_thres) out.corr.degraded|=NRTK_DEG_IONO;
        if (cell.extrap) out.corr.degraded|=NRTK_DEG_EXTRAP;
        /* too few corrected satellites (network start, outage): pass the
           nearest station through instead; back to the VRS with hysteresis */
        fallback_=out.corr.nok<(fallback_?VRS_RESUME_SATS:VRS_MIN_SATS);
        if (!fallback_) {
            out.n=nrtk_vrs_generate(vrs_,mep,&out.corr,ep.nav.get(),out.obs,MAXOBS,&out.info);
            out.degraded=out.corr.degraded;
            out.mode=1;
            return out.mode;
        }
    }
    int mode=generate_nearest(ep,out);
    if (mode&&cell.valid) out.degraded|=out.corr.degraded|NRTK_DEG_FEWSAT;
    return mode;
}

int VrsStream::generate_nearest(const NetEpoch &ep, VrsOutput &out)
{
    /* nearest station with data, healthy ones first */
    double d;
    out.n=0; out.mode=0;
    int near=nrtk_cellmgr_nearest(cm_,pos_,1,&d);
    if (near<0||!ep.sta[near]) near=nrtk_cellmgr_nearest(cm_,pos_,0,&d);
    if (near>=0&&near!=exclude_&&ep.sta[near]) {
        out.n=nrtk_vrs_passthrough(vrs_,ep.sta[near].get(),out.obs,MAXOBS,&out.info);
        out.mode=2;
        out.degraded=NRTK_DEG_NEAREST;
        out.cell.nearest=near;
    }
    return out.mode;
}

} // namespace nrtk
