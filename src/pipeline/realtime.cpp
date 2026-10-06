/*------------------------------------------------------------------------------
 * realtime.cpp : real-time core - station inputs, epoch synchronization,
 *                processing thread, VRS sessions, status and metrics
 *----------------------------------------------------------------------------*/
#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>
#include "nrtk/realtime.hpp"

namespace nrtk {

using steady=std::chrono::steady_clock;

static double secs(steady::duration d)
{
    return std::chrono::duration<double>(d).count();
}
static long long tkey(gtime_t t)
{
    return (long long)t.time*1000LL+(long long)llround(t.sec*1000.0);
}
static std::string jstr(const std::string &s)
{
    std::string o="\"";
    for (char c:s) {
        if (c=='"'||c=='\\') {o+='\\'; o+=c;}
        else if ((unsigned char)c<0x20) o+=' ';
        else o+=c;
    }
    return o+"\"";
}

/* clock ----------------------------------------------------------------------------*/
void Clock::start_replay(gtime_t t0, double speed)
{
    replay_=true; t0_=t0; speed_=speed>0.0?speed:1.0; w0_=steady::now();
}
gtime_t Clock::now() const
{
    if (!replay_) return utc2gpst(timeget());
    return timeadd(t0_,secs(steady::now()-w0_)*speed_);
}

/* epoch synchronization -------------------------------------------------------------*/
EpochSync::EpochSync(int nsta, double window_s, const Clock &clk)
    : nsta_(nsta), window_(window_s), clk_(clk), alive_(nsta)
{
}
void EpochSync::push(int sta, std::vector<obsd_t> obs)
{
    if (obs.empty()||sta<0||sta>=nsta_) return;
    gtime_t t=obs[0].time;
    long long k=tkey(t);
    double lat=timediff(clk_.now(),t);
    std::lock_guard<std::mutex> lk(mtx_);
    if (k<=last_) {late_++; return;}
    Slot &s=pend_[k];
    if (s.obs.empty()) {
        s.t=t; s.obs.resize(nsta_); s.lat.assign(nsta_,0.0); s.first=steady::now();
    }
    if (s.obs[sta].empty()) s.n++;
    s.obs[sta]=std::move(obs);
    s.lat[sta]=lat;
    alive_[sta]=steady::now();
    cv_.notify_all();
}
bool EpochSync::pop(SyncEpoch &out, int timeout_ms)
{
    auto deadline=steady::now()+std::chrono::milliseconds(timeout_ms);
    std::unique_lock<std::mutex> lk(mtx_);
    for (;;) {
        if (stop_) return false;
        auto now=steady::now();
        if (!pend_.empty()) {
            auto it=pend_.begin();
            int expected=0;
            for (auto &a:alive_) if (a.time_since_epoch().count()&&secs(now-a)<10.0) expected++;
            auto due=it->second.first+std::chrono::duration_cast<steady::duration>(
                         std::chrono::duration<double>(window_));
            if (it->second.n>=std::max(expected,1)||now>=due) {
                Slot &s=it->second;
                out.t=s.t; out.raw=std::move(s.obs); out.latency=s.lat; out.nsta=s.n;
                last_=it->first;
                pend_.erase(it);
                return true;
            }
            if (now>=deadline) return false;
            cv_.wait_until(lk,std::min(due,deadline));
            continue;
        }
        if (now>=deadline) return false;
        cv_.wait_until(lk,deadline);
    }
}
void EpochSync::shutdown()
{
    {std::lock_guard<std::mutex> lk(mtx_); stop_=true;}
    cv_.notify_all();
}

/* real-time core ----------------------------------------------------------------------*/
RealtimeCore::RealtimeCore(const nrtk_config_t &cfg, const nrtk_atx_t *atx)
    : cfg_(cfg), atx_(atx)
{
    memset(&nav_,0,sizeof(nav_));
    input_state_.assign(cfg_.nsta,"idle");
}
RealtimeCore::~RealtimeCore()
{
    stop();
    nrtk_nav_free(&nav_);
}
static bool is_replay(const char *src) {return !strncmp(src,"rinex://",8);}

bool RealtimeCore::start()
{
    net_.reset(new Network(cfg_,atx_));
    bool replay=false;
    for (int i=0;i<cfg_.nsta;i++) if (cfg_.sta[i].enabled&&is_replay(cfg_.sta[i].source)) replay=true;
    if (cfg_.nav[0]) {          /* initial ephemerides (required for replay) */
        nrtk_nav_read(cfg_.nav,&nav_);
        navdirty_=true;
    }
    sync_.reset(new EpochSync(cfg_.nsta,cfg_.sync_window_s,clk_));
    stop_=false;
    if (replay) {
        th_.emplace_back(&RealtimeCore::replay_thread,this);
    }
    else {
        for (int i=0;i<cfg_.nsta;i++) {
            if (!cfg_.sta[i].enabled||!cfg_.sta[i].source[0]) continue;
            th_.emplace_back(&RealtimeCore::input_thread,this,i);
        }
    }
    th_.emplace_back(&RealtimeCore::proc_thread,this);
    nrtk_log(NRTK_LOG_INFO,"core","started (%s, %d stations)",replay?"replay":"live",cfg_.nsta);
    return true;
}
void RealtimeCore::stop()
{
    if (stop_.exchange(true)) return;
    if (sync_) sync_->shutdown();
    bus_.shutdown();
    for (auto &t:th_) if (t.joinable()) t.join();
    th_.clear();
}

/* live input: RTKLIB stream + RTCM 3 decoder ----------------------------------*/
namespace {
struct InputCtx { RealtimeCore *core; EpochSync *sync; int ista; };
}
void RealtimeCore::add_eph(const eph_t *eph)
{
    std::lock_guard<std::mutex> lk(navmtx_);
    if (nrtk_nav_add_eph(&nav_,eph)>0) navdirty_=true;
}
void RealtimeCore::input_thread(int ista)
{
    const nrtk_station_t &sta=cfg_.sta[ista];
    /* nrtk_strsrc_t holds an rtcm_t (~0.5 MB): keep it off the thread stack */
    std::unique_ptr<nrtk_strsrc_t> srcp(new nrtk_strsrc_t);
    nrtk_strsrc_t &src=*srcp;
    while (!stop_) {
        if (!nrtk_strsrc_open(&src,sta.source)) {
            input_state_[ista]="open error";
            for (int i=0;i<50&&!stop_;i++) std::this_thread::sleep_for(std::chrono::milliseconds(100));
            continue;
        }
        input_state_[ista]="connected";
        nrtk_log(NRTK_LOG_INFO,"input","%s: stream opened",sta.id);
        auto last=steady::now();
        InputCtx ctx{this,sync_.get(),ista};
        while (!stop_) {
            int n=nrtk_strsrc_poll(&src,
                [](void *c, const obsd_t *obs, int n){
                    InputCtx *x=(InputCtx *)c;
                    std::vector<obsd_t> v(obs,obs+n);
                    for (auto &o:v) o.rcv=(uint8_t)(x->ista+1);
                    x->sync->push(x->ista,std::move(v));
                },
                [](void *c, const eph_t *eph){ ((InputCtx *)c)->core->add_eph(eph); },
                nullptr,&ctx);
            if (n>0) last=steady::now();
            else std::this_thread::sleep_for(std::chrono::milliseconds(10));
            if (n<0||secs(steady::now()-last)>30.0) {  /* reconnect */
                nrtk_log(NRTK_LOG_WARN,"input","%s: no data, reconnecting (%s)",sta.id,src.msg);
                break;
            }
        }
        nrtk_strsrc_close(&src);
        input_state_[ista]="reconnecting";
    }
    input_state_[ista]="stopped";
}

/* replay: all stations from rinex:// sources paced by the wall clock ----------*/
void RealtimeCore::replay_thread()
{
    nrtk_config_t *c=new nrtk_config_t(cfg_);
    for (int i=0;i<c->nsta;i++) {
        if (is_replay(c->sta[i].source)) snprintf(c->sta[i].obsfile,sizeof(c->sta[i].obsfile),
                                                  "%s",c->sta[i].source+8);
        else c->sta[i].obsfile[0]='\0';
    }
    c->ts=cfg_.ts; c->te=cfg_.te;
    std::unique_ptr<OfflineReader> rdp(new OfflineReader(*c));   /* holds nav_t: heap */
    OfflineReader &rd=*rdp;
    if (!rd.ok()) {
        nrtk_log(NRTK_LOG_ERROR,"replay","cannot open replay sources");
        delete c;
        return;
    }
    gtime_t t;
    std::vector<std::vector<obsd_t>> raw;
    bool first=true;
    while (!stop_&&rd.next(t,raw)) {
        if (first) {clk_.start_replay(t,cfg_.replay_speed); first=false;}
        /* wait until the epoch is due on the replay clock (+small transfer delay) */
        while (!stop_&&timediff(clk_.now(),t)<0.2) std::this_thread::sleep_for(std::chrono::milliseconds(5));
        for (int i=0;i<c->nsta;i++) {
            if (raw[i].empty()) continue;
            for (auto &o:raw[i]) o.rcv=(uint8_t)(i+1);
            for (auto &o:raw[i]) o.time=t;
            sync_->push(i,std::move(raw[i]));
            input_state_[i]="replay";
        }
    }
    nrtk_log(NRTK_LOG_INFO,"replay","end of replay data");
    rdp.reset();
    delete c;
}

std::shared_ptr<const nav_t> RealtimeCore::nav_snapshot()
{
    std::lock_guard<std::mutex> lk(navmtx_);
    if (navdirty_||!navsnap_) {
        nav_t *n=new nav_t;
        memset(n,0,sizeof(*n));
        nrtk_nav_copy(n,&nav_);
        navsnap_=std::shared_ptr<const nav_t>(n,[](const nav_t *p){
            nav_t *q=const_cast<nav_t *>(p); nrtk_nav_free(q); delete q;});
        navdirty_=false;
    }
    return navsnap_;
}

void RealtimeCore::proc_thread()
{
    SyncEpoch se;
    while (!stop_) {
        if (!sync_->pop(se,200)) continue;
        auto nav=nav_snapshot();
        if (nav->n<=0) continue;
        auto ep=net_->process(se.t,se.raw,nav,&se.latency);
        bus_.publish(ep);
        std::lock_guard<std::mutex> lk(stmtx_);
        st_.epochs++;
        st_.proc_ms_last=ep->proc_ms;
        st_.proc_ms_max=std::max(st_.proc_ms_max,ep->proc_ms);
        st_.proc_ms_sum+=ep->proc_ms;
        st_.last_t=se.t;
        for (int i=0;i<(int)se.latency.size();i++) if (!se.raw[i].empty()) {
            lat_hist_.push_back(se.latency[i]);
        }
        while (lat_hist_.size()>3000) lat_hist_.pop_front();
    }
}

RealtimeCore::Stats RealtimeCore::stats()
{
    std::lock_guard<std::mutex> lk(stmtx_);
    Stats s=st_;
    s.dropped_late=sync_?sync_->dropped_late():0;
    if (!lat_hist_.empty()) {
        std::vector<double> v(lat_hist_.begin(),lat_hist_.end());
        std::sort(v.begin(),v.end());
        s.latency_p95=v[(size_t)(0.95*(v.size()-1))];
    }
    return s;
}
void RealtimeCore::count_session_event(const char *ev)
{
    std::lock_guard<std::mutex> lk(stmtx_);
    sess_ev_[ev]++;
}

/* /status ----------------------------------------------------------------------------*/
std::string RealtimeCore::status_json(const std::vector<NtripServer::ClientStatus> &clients)
{
    std::ostringstream o;
    Stats s=stats();
    char ts[64]="";
    if (s.last_t.time) time2str(s.last_t,ts,1);
    o.setf(std::ios::fixed);
    o.precision(3);
    o<<"{\"version\":\""<<NRTK_VERSION<<"\",\"network\":"<<jstr(cfg_.name)
     <<",\"mode\":\""<<(clk_.replay()?"replay":"live")<<"\",\"last_epoch\":\""<<ts<<"\""
     <<",\"clock_tow\":"<<[this]{int w; return time2gpst(clk_.now(),&w);}()
     <<",\"epochs\":"<<s.epochs<<",\"dropped_late\":"<<s.dropped_late
     <<",\"proc_ms\":{\"last\":"<<s.proc_ms_last<<",\"max\":"<<s.proc_ms_max
     <<",\"mean\":"<<(s.epochs?s.proc_ms_sum/s.epochs:0.0)<<"}"
     <<",\"latency_p95_s\":"<<s.latency_p95<<",\"stations\":[";
    auto sts=net_->station_status();
    for (size_t i=0;i<sts.size();i++) {
        o<<(i?",":"")<<"{\"id\":"<<jstr(sts[i].id)<<",\"input\":"<<jstr(input_state_[i])
         <<",\"data\":"<<(sts[i].data?"true":"false")<<",\"healthy\":"<<(sts[i].healthy?"true":"false")
         <<",\"latency_s\":"<<sts[i].latency<<",\"age_s\":"<<std::min(sts[i].last_data_age,1E6)
         <<",\"fixrate_60s\":"<<sts[i].fixrate<<"}";
    }
    o<<"],\"baselines\":[";
    auto bls=net_->baseline_status();
    for (size_t i=0;i<bls.size();i++) {
        const auto &b=bls[i];
        o<<(i?",":"")<<"{\"pair\":\""<<cfg_.sta[b.a].id<<"-"<<cfg_.sta[b.b].id<<"\",\"km\":"<<b.length_km
         <<",\"fixrate_60s\":"<<b.fixrate_60s<<",\"fixrate_total\":"
         <<(b.nepoch?(double)b.nfixep/b.nepoch:0.0)<<",\"ttff_s\":"<<b.ttff<<",\"slips\":"<<b.nslip
         <<",\"refix\":"<<b.nrefix<<",\"iono_rate_m_per_min\":"<<b.iono_rate<<"}";
    }
    o<<"],\"sessions\":[";
    for (size_t i=0;i<clients.size();i++) {
        const auto &c=clients[i];
        o<<(i?",":"")<<"{\"id\":"<<c.info.id<<",\"mount\":"<<jstr(c.info.mount)<<",\"user\":"
         <<jstr(c.info.user)<<",\"peer\":"<<jstr(c.info.peer)<<",\"connected_s\":"<<c.connected_s
         <<",\"bytes\":"<<c.bytes_sent<<",\"gga\":"<<c.ngga<<",\"session\":"
         <<(c.session_json.empty()?"null":c.session_json)<<"}";
    }
    o<<"],\"session_events\":{";
    {
        std::lock_guard<std::mutex> lk(stmtx_);
        int k=0;
        for (auto &kv:sess_ev_) o<<(k++?",":"")<<jstr(kv.first)<<":"<<kv.second;
    }
    o<<"}}";
    return o.str();
}

/* /metrics (Prometheus text format) ---------------------------------------------------*/
std::string RealtimeCore::metrics(const std::vector<NtripServer::ClientStatus> &clients)
{
    std::ostringstream o;
    Stats s=stats();
    o.setf(std::ios::fixed); o.precision(4);
    o<<"# HELP nrtk_epochs_total Processed network epochs\n# TYPE nrtk_epochs_total counter\n"
     <<"nrtk_epochs_total "<<s.epochs<<"\n";
    o<<"# HELP nrtk_epoch_processing_ms Processing time of the last epoch\n# TYPE nrtk_epoch_processing_ms gauge\n"
     <<"nrtk_epoch_processing_ms "<<s.proc_ms_last<<"\n";
    o<<"# HELP nrtk_station_latency_p95_seconds 95th percentile data latency\n# TYPE nrtk_station_latency_p95_seconds gauge\n"
     <<"nrtk_station_latency_p95_seconds "<<s.latency_p95<<"\n";
    o<<"# HELP nrtk_epochs_late_total Station epochs arriving after the sync window\n# TYPE nrtk_epochs_late_total counter\n"
     <<"nrtk_epochs_late_total "<<s.dropped_late<<"\n";
    o<<"# HELP nrtk_station_healthy Station health (1 healthy)\n# TYPE nrtk_station_healthy gauge\n";
    auto sts=net_->station_status();
    for (auto &x:sts) o<<"nrtk_station_healthy{station=\""<<x.id<<"\"} "<<(x.healthy?1:0)<<"\n";
    o<<"# HELP nrtk_station_latency_seconds Data latency of the last epoch\n# TYPE nrtk_station_latency_seconds gauge\n";
    for (auto &x:sts) o<<"nrtk_station_latency_seconds{station=\""<<x.id<<"\"} "<<x.latency<<"\n";
    o<<"# HELP nrtk_baseline_fixrate Fixed epochs over the last 60 s\n# TYPE nrtk_baseline_fixrate gauge\n";
    auto bls=net_->baseline_status();
    for (auto &b:bls) o<<"nrtk_baseline_fixrate{baseline=\""<<cfg_.sta[b.a].id<<"-"<<cfg_.sta[b.b].id<<"\"} "
                      <<b.fixrate_60s<<"\n";
    o<<"# HELP nrtk_baseline_iono_rate DD ionosphere rate (m/min)\n# TYPE nrtk_baseline_iono_rate gauge\n";
    for (auto &b:bls) o<<"nrtk_baseline_iono_rate{baseline=\""<<cfg_.sta[b.a].id<<"-"<<cfg_.sta[b.b].id<<"\"} "
                      <<b.iono_rate<<"\n";
    o<<"# HELP nrtk_sessions Connected NTRIP clients\n# TYPE nrtk_sessions gauge\n"
     <<"nrtk_sessions "<<clients.size()<<"\n";
    o<<"# HELP nrtk_session_events_total Session events\n# TYPE nrtk_session_events_total counter\n";
    {
        std::lock_guard<std::mutex> lk(stmtx_);
        for (auto &kv:sess_ev_) o<<"nrtk_session_events_total{event=\""<<kv.first<<"\"} "<<kv.second<<"\n";
    }
    return o.str();
}

/* VRS session ---------------------------------------------------------------------------*/
VrsSession::VrsSession(RealtimeCore &core, const ClientInfo &ci, int msm, bool nearest_only)
    : core_(core), ci_(ci), msm_(msm), nearest_only_(nearest_only)
{
    created_=last_ok_=steady::now();
    staid_base_=(int)(1000+(ci.id*37)%3000);
    core_.count_session_event("connect");
    nrtk_log(NRTK_LOG_INFO,"session","id=%llu mount=%s user=%s peer=%s created",
             (unsigned long long)ci.id,ci.mount.c_str(),ci.user.c_str(),ci.peer.c_str());
}
VrsSession::~VrsSession()
{
    if (enc_init_) nrtk_rtcmout_free(&enc_);
    core_.count_session_event("close");
    nrtk_log(NRTK_LOG_INFO,"session","id=%llu closed after %.0f s, %ld epochs sent, %ld rebuilds",
             (unsigned long long)ci_.id,secs(steady::now()-created_),nout_,nrebuild_);
}
void VrsSession::build_vrs(const double *pos)
{
    double p[3];
    for (int i=0;i<3;i++) p[i]=floor(pos[i]+0.5);       /* 1 m grid */
    nrtk_cell_t cell;
    nrtk_cellmgr_select(core_.network().cellmgr(),p,&cell);
    memcpy(vpos_,p,sizeof(vpos_));
    vrs_.reset(new VrsStream(core_.network(),vpos_,-1));
    for (int k=0;k<3;k++) tri_[k]=cell.tri[k];
    staid_=have_vrs_?((staid_+1)&0xFFF):staid_base_;
    if (!enc_init_) {nrtk_rtcmout_init(&enc_,staid_,msm_,vpos_,1); enc_init_=true;}
    else nrtk_rtcmout_restart(&enc_,staid_,vpos_);
    enc_.non_physical=1;
    enc_.sysmask=core_.config().vrs.sysmask;
    cur_station_=-1;
    if (have_vrs_) {nrebuild_++; core_.count_session_event("rebuild");}
    else core_.count_session_event("vrs_create");
    have_vrs_=true;
    double llh[3]; ecef2pos(vpos_,llh);
    nrtk_log(NRTK_LOG_INFO,"session","id=%llu VRS %s at %.6f %.6f %.1f staid=%d master=%s%s",
             (unsigned long long)ci_.id,nrebuild_?"rebuilt":"created",llh[0]*R2D,llh[1]*R2D,llh[2],
             staid_,cell.valid?core_.config().sta[cell.master].id:"-",cell.extrap?" (extrapolated)":"");
}
void VrsSession::on_gga(const GgaFix &fix)
{
    if (!fix.valid) return;
    fix_=fix; have_fix_=true;
    if (!have_vrs_) {build_vrs(fix.ecef); return;}
    /* rebuild: rover moved more than vrs_rebuild_km or left the triangle */
    double d=nrtk_dist(fix.ecef,vpos_)/1E3;
    bool moved=d>core_.config().vrs_rebuild_km;
    bool cellchg=false;
    if (!moved&&!nearest_only_) {
        nrtk_cell_t cell;
        nrtk_cellmgr_select(core_.network().cellmgr(),fix.ecef,&cell);
        cellchg=cell.valid&&tri_[0]>=0&&
                !(std::is_permutation(cell.tri,cell.tri+3,tri_));
    }
    if (moved||cellchg) build_vrs(fix.ecef);
}
bool VrsSession::on_epoch(const void *epoch, std::vector<uint8_t> &out)
{
    const NetEpoch &ep=*(const NetEpoch *)epoch;
    const nrtk_config_t &cfg=core_.config();
    out.clear();
    nep_++;
    if (!have_vrs_) return true;                 /* waiting for a valid GGA */

    VrsOutput *vo=new VrsOutput;
    int mode=nearest_only_?vrs_->generate_nearest(ep,*vo):vrs_->generate(ep,*vo);
    if (mode==0||vo->n<=0) {
        delete vo;
        if (secs(steady::now()-last_ok_)>cfg.nodata_close_s) {
            nrtk_log(NRTK_LOG_WARN,"session","id=%llu no VRS for %.0f s: closing",
                     (unsigned long long)ci_.id,cfg.nodata_close_s);
            core_.count_session_event("close_nodata");
            return false;
        }
        return true;
    }
    /* reference station identity: VRS, or the passed-through station */
    int station=mode==1?-1:vo->cell.nearest;
    bool restart=station!=cur_station_;
    if (restart) {
        if (station<0) nrtk_rtcmout_restart(&enc_,staid_,vpos_);
        else {
            double pos[3];
            nrtk_sta_pos(cfg.sta+station,ep.t,pos);
            nrtk_rtcmout_restart(&enc_,cfg.sta[station].rtcm_id?cfg.sta[station].rtcm_id:staid_,pos);
        }
        enc_.non_physical=station<0;
        if (cur_station_!=-1||station!=-1) {
            core_.count_session_event(station<0?"vrs_resume":"degrade_nearest");
            nrtk_log(NRTK_LOG_INFO,"session","id=%llu %s",(unsigned long long)ci_.id,
                     station<0?"VRS resumed":"degraded: nearest station pass-through");
        }
        cur_station_=station;
    }
    uint8_t *buf=new uint8_t[65536];
    int n=0;
    if (restart||nout_%10==0) n+=nrtk_rtcmout_station(&enc_,buf+n,65536-n);
    n+=nrtk_rtcmout_obs(&enc_,ep.t,vo->obs,vo->n,buf+n,65536-n);
    if (cfg.eph_out&&ep.nav) {                  /* one satellite ephemeris per epoch */
        int sat=vo->obs[eph_cursor_++%vo->n].sat;
        n+=nrtk_rtcmout_eph(&enc_,ep.nav.get(),sat,buf+n,65536-n);
    }
    out.assign(buf,buf+n);
    delete[] buf;
    last_ok_=steady::now();
    nout_++;
    last_mode_=mode; last_nsat_=vo->n; last_deg_=vo->degraded; last_t_=ep.t;
    last_age_=timediff(core_.clock().now(),ep.t);
    last_master_=mode==1?cfg.sta[vo->cell.master].id:cfg.sta[vo->cell.nearest].id;
    delete vo;
    return true;
}
std::string VrsSession::status_json() const
{
    std::ostringstream o;
    double llh[3]={0,0,0};
    if (have_vrs_) ecef2pos(vpos_,llh);
    o.setf(std::ios::fixed); o.precision(7);
    o<<"{\"vrs\":"<<(have_vrs_?"true":"false")<<",\"lat\":"<<llh[0]*R2D<<",\"lon\":"<<llh[1]*R2D;
    o.precision(2);
    o<<",\"h\":"<<llh[2]<<",\"staid\":"<<staid_<<",\"mode\":\""
     <<(last_mode_==1?"vrs":(last_mode_==2?"nearest":"none"))<<"\",\"master\":"<<jstr(last_master_)
     <<",\"nsat\":"<<last_nsat_<<",\"degraded\":"<<last_deg_<<",\"age_s\":"<<last_age_
     <<",\"epochs_sent\":"<<nout_<<",\"rebuilds\":"<<nrebuild_<<"}";
    return o.str();
}

} // namespace nrtk
