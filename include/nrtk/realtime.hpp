/*------------------------------------------------------------------------------
 * nrtk/realtime.hpp : real-time core (spec ch.11)
 *
 *   station inputs (thread per station: NTRIP/TCP/serial/file via RTKLIB
 *   stream + RTCM 3 decoder, or rinex:// replay paced by the wall clock)
 *     -> EpochSync (collects one GPS epoch from all stations, waits at most
 *        sync_window_s after the first arrival; late stations are absent)
 *     -> processing thread: Network::process (baselines on a thread pool)
 *     -> EpochBus (latest NetEpoch + 30 s history = corrections cache)
 *     -> one VrsSession per NTRIP client (GGA -> VRS -> RTCM 3)
 *----------------------------------------------------------------------------*/
#ifndef NRTK_REALTIME_HPP
#define NRTK_REALTIME_HPP

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include "nrtk/pipeline.hpp"
#include "nrtk/ntrip_front.hpp"
#include "nrtk/rtcm_out.h"

namespace nrtk {

/* time base: real GPST or a replay clock running from the data start --------*/
class Clock {
public:
    void start_replay(gtime_t t0, double speed);
    bool replay() const {return replay_;}
    gtime_t now() const;                 /* GPST */
    double  speed() const {return speed_;}
private:
    bool replay_=false;
    gtime_t t0_{};
    double speed_=1.0;
    std::chrono::steady_clock::time_point w0_;
};

struct SyncEpoch {
    gtime_t t{};
    std::vector<std::vector<obsd_t>> raw;    /* per station, empty if absent */
    std::vector<double> latency;             /* arrival - epoch time (s) */
    int nsta=0;
};

class EpochSync {
public:
    EpochSync(int nsta, double window_s, const Clock &clk);
    void push(int sta, std::vector<obsd_t> obs);
    bool pop(SyncEpoch &out, int timeout_ms);
    void shutdown();
    long dropped_late() const {return late_;}
private:
    struct Slot {
        gtime_t t{};
        std::vector<std::vector<obsd_t>> obs;
        std::vector<double> lat;
        int n=0;
        std::chrono::steady_clock::time_point first;
    };
    int nsta_;
    double window_;
    const Clock &clk_;
    std::map<long long,Slot> pend_;          /* key: ms of GPST */
    std::vector<std::chrono::steady_clock::time_point> alive_;
    long long last_=-1;
    std::atomic<long> late_{0};
    std::mutex mtx_;
    std::condition_variable cv_;
    bool stop_=false;
};

class RealtimeCore;

/* VRS session (one NTRIP connection) ------------------------------------------*/
class VrsSession : public StreamSession {
public:
    VrsSession(RealtimeCore &core, const ClientInfo &ci, int msm, bool nearest_only);
    ~VrsSession() override;
    void on_gga(const GgaFix &fix) override;
    bool on_epoch(const void *epoch, std::vector<uint8_t> &out) override;
    std::string status_json() const override;
private:
    void build_vrs(const double *pos);
    RealtimeCore &core_;
    ClientInfo ci_;
    int msm_;
    bool nearest_only_;
    GgaFix fix_;
    bool have_fix_=false;
    double vpos_[3]={0,0,0};
    bool have_vrs_=false;
    int staid_=0,staid_base_=0,cur_station_=-1,tri_[3]={-1,-1,-1};
    std::unique_ptr<VrsStream> vrs_;
    nrtk_rtcmout_t enc_{};
    bool enc_init_=false;
    std::chrono::steady_clock::time_point last_ok_,created_;
    long nep_=0,nout_=0,nrebuild_=0;
    int last_mode_=0,last_nsat_=0,last_deg_=0;
    double last_age_=0.0;
    std::string last_master_;
    gtime_t last_t_{};
    int eph_cursor_=0;
};

/* real-time core ----------------------------------------------------------------*/
class RealtimeCore {
public:
    RealtimeCore(const nrtk_config_t &cfg, const nrtk_atx_t *atx);
    ~RealtimeCore();
    bool start();
    void stop();
    const nrtk_config_t &config() const {return cfg_;}
    Network &network() {return *net_;}
    const Clock &clock() const {return clk_;}
    EpochBus<NetEpoch> &bus() {return bus_;}
    std::shared_ptr<const nav_t> nav_snapshot();
    std::string status_json(const std::vector<NtripServer::ClientStatus> &clients);
    std::string metrics(const std::vector<NtripServer::ClientStatus> &clients);
    void count_session_event(const char *ev);
    struct Stats {
        long epochs=0,dropped_late=0;
        double proc_ms_last=0.0,proc_ms_max=0.0,proc_ms_sum=0.0;
        double latency_p95=0.0;
        gtime_t last_t{};
    };
    Stats stats();
private:
    void input_thread(int ista);
    void replay_thread();
    void proc_thread();
    void add_eph(const eph_t *eph);
    nrtk_config_t cfg_;
    const nrtk_atx_t *atx_;
    std::unique_ptr<Network> net_;
    Clock clk_;
    std::unique_ptr<EpochSync> sync_;
    EpochBus<NetEpoch> bus_;
    std::vector<std::thread> th_;
    std::atomic<bool> stop_{false};
    std::mutex navmtx_;
    nav_t nav_{};
    std::shared_ptr<const nav_t> navsnap_;
    bool navdirty_=true;
    std::mutex stmtx_;
    Stats st_;
    std::deque<double> lat_hist_;
    std::map<std::string,long> sess_ev_;
    std::vector<std::string> input_state_;
};

} // namespace nrtk

#endif
