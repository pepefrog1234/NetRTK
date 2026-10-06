/*------------------------------------------------------------------------------
 * nrtk/pipeline.hpp : processing pipeline shared by nrtk-gen (offline) and
 *                     nrtk-svr (real time)
 *
 *   raw obs per station --prep--> prepared epochs --netpos--> baseline epochs
 *   (all station pairs <= max_baseline_km, thread pool)       |
 *                                                             v
 *   VRS request (position) --cellmgr--> cell --atmos--> corrections --vrsgen-->
 *   VRS observations --rtcm_out--> RTCM 3 / RINEX
 *----------------------------------------------------------------------------*/
#ifndef NRTK_PIPELINE_HPP
#define NRTK_PIPELINE_HPP

#include <atomic>
#include <condition_variable>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include "nrtk/config.h"
#include "nrtk/antenna.h"
#include "nrtk/prep.h"
#include "nrtk/netpos.h"
#include "nrtk/atmos.h"
#include "nrtk/cellmgr.h"
#include "nrtk/vrsgen.h"
#include "nrtk/ingest.h"

namespace nrtk {

/* simple fixed-size thread pool ----------------------------------------------*/
class ThreadPool {
public:
    explicit ThreadPool(int n);
    ~ThreadPool();
    /* run f(i) for i in [0,n) and wait */
    void parallel_for(int n, const std::function<void(int)> &f);
private:
    void worker();
    std::vector<std::thread> th_;
    std::mutex mtx_;
    std::condition_variable cv_,done_cv_;
    std::function<void(int)> job_;
    int next_=0,total_=0,pending_=0;
    long gen_=0;
    bool stop_=false;
};

/* results of one network epoch (immutable once published) -------------------*/
struct NetEpoch {
    gtime_t t{};
    std::vector<std::shared_ptr<nrtk_staepoch_t>> sta;   /* prepared, null: absent */
    std::map<std::pair<int,int>,nrtk_blepoch_t> bl;      /* key (master side, aux side) */
    std::shared_ptr<const nav_t> nav;                    /* ephemeris snapshot */
    std::vector<int> healthy;                            /* station health */
    std::vector<double> latency;                         /* data latency (s) */
    double proc_ms=0.0;                                  /* processing time */
    const nrtk_blepoch_t *find(int a, int b, int *flip) const;
};

struct StationStatus {
    std::string id;
    bool   data=false, healthy=false, recovering=false;
    double latency=0.0, fixrate=0.0;
    double last_data_age=1E9;
};

struct BaselineStatus {
    int    a=0,b=0;
    double length_km=0.0, fixrate_60s=0.0, ratio=0.0, ttff=-1.0;
    long   nepoch=0,nfixep=0,nslip=0,nrefix=0;
    int    nsat=0,nfixed=0;
    double iono_rate=0.0;
};

/* network processing ------------------------------------------------------------*/
class Network {
public:
    Network(const nrtk_config_t &cfg, const nrtk_atx_t *atx);
    ~Network();
    Network(const Network &)=delete;
    Network &operator=(const Network &)=delete;

    /* process one synchronized epoch. raw[i]: observations of station i
     * (empty: absent). latency[i]: data latency (real time, else 0). */
    std::shared_ptr<NetEpoch> process(gtime_t t,
                                      const std::vector<std::vector<obsd_t>> &raw,
                                      std::shared_ptr<const nav_t> nav,
                                      const std::vector<double> *latency=nullptr);
    const nrtk_config_t &config() const {return cfg_;}
    nrtk_cellmgr_t *cellmgr() {return cm_;}
    const nrtk_ionomon_t *ionomon() const {return mon_;}
    std::vector<StationStatus> station_status() const;
    std::vector<BaselineStatus> baseline_status() const;
    int nsta() const {return cfg_.nsta;}

    long closure_flags() const {return nclosure_;}
private:
    void update_health(NetEpoch &ep);
    void closure_check(NetEpoch &ep);
    long nclosure_=0;
    nrtk_config_t cfg_;
    const nrtk_atx_t *atx_;
    std::vector<const nrtk_ant_t *> ant_;
    struct BL { int m,a; nrtk_baseline_t *bl; };
    std::vector<BL> bls_;
    nrtk_cellmgr_t *cm_=nullptr;
    nrtk_ionomon_t *mon_=nullptr;
    std::unique_ptr<ThreadPool> pool_;
    gtime_t tstart_{};
    std::vector<StationStatus> st_;
    std::vector<gtime_t> last_data_;
    mutable std::mutex mtx_;
};

/* one VRS stream (offline job or NTRIP session) --------------------------------*/
struct VrsOutput {
    int    n=0;
    obsd_t obs[MAXOBS];
    nrtk_vcorr_t corr{};
    nrtk_vrs_info_t info{};
    nrtk_cell_t cell{};
    int    mode=0;               /* 0: none, 1: VRS, 2: nearest pass-through */
    int    degraded=0;
    double iono_activity=0.0;
};

class VrsStream {
public:
    /* exclude: station index left out of the network (leave-one-out), -1 none */
    VrsStream(Network &net, const double *pos, int exclude=-1);
    ~VrsStream();
    VrsStream(const VrsStream &)=delete;
    VrsStream &operator=(const VrsStream &)=delete;
    /* generate one epoch. returns mode (0: nothing) */
    int generate(const NetEpoch &ep, VrsOutput &out);
    /* nearest healthy station pass-through only (NEAR mountpoint) */
    int generate_nearest(const NetEpoch &ep, VrsOutput &out);
    const double *pos() const {return pos_;}
    void reset();
private:
    Network &net_;
    double pos_[3];
    int exclude_;
    nrtk_cellmgr_t *cm_=nullptr;   /* own cell manager (leave-one-out) */
    nrtk_vrs_t *vrs_=nullptr;
    std::vector<int> health_;
    bool fallback_=false;          /* too few corrected satellites: nearest station */
    nrtk_zwd_state_t zwd_{};        /* smoothed wet anomaly of the cell master (atmos) */
};

/* offline input: epoch-aligned reading of the stations' RINEX files ----------*/
class OfflineReader {
public:
    explicit OfflineReader(const nrtk_config_t &cfg);
    ~OfflineReader();
    bool ok() const {return ok_;}
    const nav_t &nav() const {return nav_;}
    /* shared handle to the (static) navigation data of the whole run */
    std::shared_ptr<const nav_t> nav_ptr() const {return std::shared_ptr<const nav_t>(&nav_,[](const nav_t *){});}
    /* next epoch on the time grid; false at the end */
    bool next(gtime_t &t, std::vector<std::vector<obsd_t>> &raw);
    double interval() const {return tint_;}
private:
    const nrtk_config_t &cfg_;
    std::vector<nrtk_rnxsrc_t> src_;
    std::vector<int> open_;
    nav_t nav_{};
    gtime_t t_{},te_{};
    double tint_=1.0;
    bool ok_=false,first_=true;
};

/* VRS writers for offline jobs (RINEX 3.04, RTCM 3, diagnostics CSV) --------*/
class VrsWriter {
public:
    VrsWriter(const nrtk_config_t &cfg, const std::string &dir, const std::string &name,
              const double *pos, int staid, const nav_t *nav);
    ~VrsWriter();
    void write(const NetEpoch &ep, const VrsOutput &out);
    void close();
private:
    const nrtk_config_t &cfg_;
    FILE *rnx_=nullptr,*rtcm_=nullptr,*diag_=nullptr,*epo_=nullptr;
    rnxopt_t ropt_;
    struct RtcmState;
    RtcmState *rs_=nullptr;
    double pos_[3];
    long nep_=0;
};

/* time helpers */
std::string time_str(gtime_t t, int ndec=1);
int  mkdirs(const std::string &dir);

} // namespace nrtk

#endif
