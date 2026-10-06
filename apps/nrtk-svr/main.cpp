/*------------------------------------------------------------------------------
 * nrtk-svr : real-time network RTK server (milestones M5/M6)
 *
 *   nrtk-svr -c config.yaml [-p port] [-v level] [--json-log]
 *
 * Mountpoints: <mount_msm7> (VRS, MSM7), <mount_msm4> (VRS, MSM4) and,
 * if enabled, <mount_nearest> (nearest station pass-through). HTTP GET
 * /status returns JSON, /metrics Prometheus text.
 *----------------------------------------------------------------------------*/
#include <csignal>
#include <cstdio>
#include <cstring>
#include <memory>
#include "nrtk/realtime.hpp"
#include "../common_cli.hpp"

using namespace nrtk;

static volatile sig_atomic_t g_stop=0;
static void on_signal(int) {g_stop=1;}

static const char *usage=
"usage: nrtk-svr -c config.yaml [-p port] [-v level] [--json-log]\n";

int main(int argc, char **argv)
{
    std::string conf;
    int port=-1,verbose=NRTK_LOG_INFO,json=0;
    for (int i=1;i<argc;i++) {
        std::string a=argv[i];
        if (a=="-c"&&i+1<argc) conf=argv[++i];
        else if (a=="-p"&&i+1<argc) port=atoi(argv[++i]);
        else if (a=="-v"&&i+1<argc) verbose=atoi(argv[++i]);
        else if (a=="--json-log") json=1;
        else {fputs(usage,stderr); return 2;}
    }
    if (conf.empty()) {fputs(usage,stderr); return 2;}
    CliOpts o; o.conf=conf; o.verbose=verbose;
    nrtk_config_t *cfg=new nrtk_config_t;
    nrtk_atx_t atx;
    if (!load_config(o,*cfg,atx)) return 2;
    if (port>=0) cfg->port=port;
    nrtk_log_json(json||cfg->log_json);
    if (cfg->log_file[0]&&!freopen(cfg->log_file,"a",stderr)) {
        fprintf(stdout,"cannot open log file %s\n",cfg->log_file);
    }
    signal(SIGINT,on_signal); signal(SIGTERM,on_signal);
    signal(SIGPIPE,SIG_IGN);

    RealtimeCore core(*cfg,&atx);
    core.start();

    /* mountpoints and sourcetable */
    double lat=0,lon=0;
    for (int i=0;i<cfg->nsta;i++) {double p[3]; ecef2pos(cfg->sta[i].pos,p); lat+=p[0]*R2D/cfg->nsta; lon+=p[1]*R2D/cfg->nsta;}
    std::string sys;
    const char *names[]={"GPS","GAL","BDS","QZS"};
    for (int i=0;i<NRTK_NSYS;i++) if (cfg->vrs.sysmask&(1<<i)) sys+=std::string(sys.empty()?"":"+")+names[i];
    std::vector<MountPoint> mounts;
    auto mk=[&](const char *name, const char *id, const char *det, bool gga){
        MountPoint m; m.name=name; m.identifier=id; m.format="RTCM 3.3"; m.format_details=det;
        m.nav_system=sys; m.network=cfg->name; m.country="TWN"; m.lat=lat; m.lon=lon; m.need_gga=gga;
        mounts.push_back(m);
    };
    mk(cfg->mount_msm7,"VRS","1006(10),1033(10),1077(1),1097(1),1127(1),1117(1)",true);
    mk(cfg->mount_msm4,"VRS","1006(10),1033(10),1074(1),1094(1),1124(1),1114(1)",true);
    if (cfg->enable_near) mk(cfg->mount_near,"NEAREST","1006(10),1033(10),1077(1),1097(1),1127(1),1117(1)",true);

    NtripServerOptions so;
    so.port=cfg->port;
    so.gga_timeout_s=cfg->gga_timeout_s;
    so.allow_anonymous=cfg->allow_anonymous;
    so.network_name=cfg->name;
    for (int i=0;i<cfg->nuser;i++) {
        UserAccount u; u.user=cfg->user[i].user; u.pass=cfg->user[i].pass; u.max_conn=cfg->user[i].max_conn;
        so.users.push_back(u);
    }
    std::string m7=cfg->mount_msm7,m4=cfg->mount_msm4,mn=cfg->mount_near;
    NtripServer *srvp=nullptr;
    auto factory=[&](const ClientInfo &ci)->std::unique_ptr<StreamSession>{
        if (ci.mount==m7) return std::unique_ptr<StreamSession>(new VrsSession(core,ci,7,false));
        if (ci.mount==m4) return std::unique_ptr<StreamSession>(new VrsSession(core,ci,4,false));
        if (ci.mount==mn&&cfg->enable_near) return std::unique_ptr<StreamSession>(new VrsSession(core,ci,7,true));
        return nullptr;
    };
    auto wait_epoch=[&](uint64_t *seq,int timeout_ms)->std::shared_ptr<const void>{
        return core.bus().wait(seq,timeout_ms);
    };
    auto http=[&](const std::string &path,std::string &ctype,std::string &body)->bool{
        auto clients=srvp?srvp->clients():std::vector<NtripServer::ClientStatus>();
        if (path=="/status") {ctype="application/json"; body=core.status_json(clients); return true;}
        if (path=="/metrics") {ctype="text/plain; version=0.0.4"; body=core.metrics(clients); return true;}
        return false;
    };
    NtripServer srv(so,mounts,factory,wait_epoch,http);
    srvp=&srv;
    if (!srv.start()) {
        nrtk_log(NRTK_LOG_ERROR,"svr","cannot start NTRIP server on port %d",cfg->port);
        core.stop();
        return 1;
    }
    nrtk_log(NRTK_LOG_INFO,"svr","NTRIP caster on port %d, mountpoints %s %s%s%s",srv.port(),
             m7.c_str(),m4.c_str(),cfg->enable_near?" ":"",cfg->enable_near?mn.c_str():"");
    int k=0;
    while (!g_stop) {
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        if (++k%300==0) {   /* every minute */
            auto st=core.stats();
            nrtk_log(NRTK_LOG_INFO,"svr","epochs %ld, proc %.1f ms (max %.1f), latency p95 %.2f s, clients %zu",
                     st.epochs,st.proc_ms_last,st.proc_ms_max,st.latency_p95,srv.clients().size());
        }
    }
    nrtk_log(NRTK_LOG_INFO,"svr","shutting down");
    srv.stop();
    core.stop();
    nrtk_atx_free(&atx);
    delete cfg;
    return 0;
}
