/*------------------------------------------------------------------------------
 * nrtk/ntrip_front.hpp : NTRIP caster front end (v1 and v2) with per-connection
 *                        sessions and rover GGA handling (spec ch.8 / ch.11)
 *
 * The server knows nothing about VRS. For every client that requests a
 * mountpoint it asks a factory for a StreamSession; the connection thread
 *   - reads NMEA GGA sentences sent upstream by the rover and passes valid
 *     fixes to StreamSession::on_gga(),
 *   - waits for new epochs from the EpochSource and writes the bytes produced
 *     by StreamSession::on_epoch() to the client,
 *   - closes the connection when on_epoch() returns false, when no valid GGA
 *     was received within gga_timeout_s for a mountpoint that needs GGA, or
 *     when the client disconnects.
 * Non-NTRIP HTTP paths (e.g. /status, /metrics) are delegated to a handler.
 *----------------------------------------------------------------------------*/
#ifndef NRTK_NTRIP_FRONT_HPP
#define NRTK_NTRIP_FRONT_HPP

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include "rtklib.h"

namespace nrtk {

/* rover position from an NMEA GGA sentence ----------------------------------*/
struct GgaFix {
    bool    valid=false;
    int     quality=0;          /* GGA fix quality (0 invalid, 1 GPS, 2 DGPS, 4 fix, 5 float...) */
    int     nsat=0;
    double  hdop=0.0;
    double  lat=0.0,lon=0.0;    /* rad */
    double  h_msl=0.0;          /* orthometric height (m) */
    double  geoid=0.0;          /* geoid separation (m), from GGA or EGM96 */
    bool    geoid_from_gga=false;
    double  h_ell=0.0;          /* ellipsoidal height = h_msl + geoid (m) */
    double  tod=0.0;            /* UTC time of day (s) */
    double  ecef[3]={0,0,0};
};

/* Parse one NMEA sentence ($GPGGA/$GNGGA/$..GGA). Checks the checksum when
 * present. Missing geoid separation is replaced by EGM96 (RTKLIB geoidh).
 * Returns true for a syntactically valid GGA (fix.valid: quality>=1). */
bool parse_gga(const std::string &sentence, GgaFix &fix);

/* epoch notifications from the processing core ------------------------------*/
template<class Epoch> class EpochBus {
public:
    /* publish a new epoch (processing thread) */
    void publish(std::shared_ptr<const Epoch> ep) {
        {std::lock_guard<std::mutex> lk(mtx_); last_=std::move(ep); seq_++;}
        cv_.notify_all();
    }
    /* wait for an epoch newer than *seq (connection threads). returns null on
     * timeout or shutdown; updates *seq */
    std::shared_ptr<const Epoch> wait(uint64_t *seq, int timeout_ms) {
        std::unique_lock<std::mutex> lk(mtx_);
        if (!cv_.wait_for(lk,std::chrono::milliseconds(timeout_ms),
                          [&]{return stop_||seq_!=*seq;})) return nullptr;
        if (stop_) return nullptr;
        *seq=seq_;
        return last_;
    }
    uint64_t seq() const {std::lock_guard<std::mutex> lk(mtx_); return seq_;}
    void shutdown() {{std::lock_guard<std::mutex> lk(mtx_); stop_=true;} cv_.notify_all();}
private:
    mutable std::mutex mtx_;
    std::condition_variable cv_;
    std::shared_ptr<const Epoch> last_;
    uint64_t seq_=0;
    bool stop_=false;
};

/* one client connection's data producer -------------------------------------*/
class StreamSession {
public:
    virtual ~StreamSession()=default;
    /* valid GGA received from the rover (wall-clock receive time in GPST) */
    virtual void on_gga(const GgaFix &fix)=0;
    /* bytes to send for a new epoch (empty: nothing to send).
     * return false to close the connection. */
    virtual bool on_epoch(const void *epoch, std::vector<uint8_t> &out)=0;
    /* session description for /status (JSON object) */
    virtual std::string status_json() const=0;
};

struct MountPoint {
    std::string name;           /* e.g. VRS_MSM7 */
    std::string identifier;     /* sourcetable identifier (city / description) */
    std::string format;         /* RTCM 3.3 */
    std::string format_details; /* 1006(10),1033(10),1077(1),1087(1),... */
    std::string nav_system;     /* GPS+GAL+BDS+QZS */
    std::string network;        /* network name */
    std::string country;        /* TWN */
    double lat=0.0,lon=0.0;     /* deg, approximate */
    bool   need_gga=true;       /* sourcetable nmea=1 */
    std::string generator="nrtk";
};

struct UserAccount {
    std::string user,pass;
    int max_conn=4;
};

struct NtripServerOptions {
    int    port=2101;
    std::string bind_addr="0.0.0.0";
    double gga_timeout_s=15.0;  /* close if no valid GGA (need_gga mounts) */
    int    max_clients=256;
    bool   allow_anonymous=true;
    std::vector<UserAccount> users;
    std::string caster_host="nrtk";
    std::string operator_name="nrtk";
    std::string network_name="NRTK";
};

/* per-connection information given to the session factory */
struct ClientInfo {
    uint64_t    id=0;
    std::string mount,user,peer,agent;
    int         ntrip_version=1;        /* 1 or 2 */
};

using SessionFactory=std::function<std::unique_ptr<StreamSession>(const ClientInfo &)>;
/* HTTP handler for non-mountpoint paths. returns false for 404 */
using HttpHandler=std::function<bool(const std::string &path, std::string &content_type,
                                     std::string &body)>;

class NtripServer {
public:
    NtripServer(const NtripServerOptions &opt, std::vector<MountPoint> mounts,
                SessionFactory factory, std::function<std::shared_ptr<const void>(uint64_t *seq,int timeout_ms)> wait_epoch,
                HttpHandler http=nullptr);
    ~NtripServer();
    bool start();               /* bind and start the accept thread */
    void stop();                /* close all connections and join threads */
    int  port() const;          /* bound port (useful with port 0 in tests) */

    struct ClientStatus {
        ClientInfo info;
        double connected_s=0.0;
        long   bytes_sent=0;
        int    ngga=0;
        std::string session_json;
    };
    std::vector<ClientStatus> clients() const;
    std::string sourcetable() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/* base64 helpers (exposed for tests) */
std::string base64_encode(const std::string &in);
std::string base64_decode(const std::string &in);

} // namespace nrtk

#endif
