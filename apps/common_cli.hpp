/* shared command line handling for nrtk tools */
#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "nrtk/config.h"
#include "nrtk/antenna.h"

namespace nrtk {

struct CliOpts {
    std::string conf, outdir, ts, te, job;
    std::string datum;          /* nrtk-gen: output datum of the RTCM 1006 position */
    double dur=0.0, tint=0.0;
    int verbose=NRTK_LOG_INFO;
    bool satcsv=false;
};

inline bool parse_cli(int argc, char **argv, CliOpts &o, const char *usage)
{
    for (int i=1;i<argc;i++) {
        std::string a=argv[i];
        auto val=[&](){ if (i+1>=argc) {fprintf(stderr,"%s needs a value\n",a.c_str()); exit(2);} return std::string(argv[++i]); };
        if (a=="-c") o.conf=val();
        else if (a=="-o") o.outdir=val();
        else if (a=="-ts") o.ts=val();
        else if (a=="-te") o.te=val();
        else if (a=="-d") o.dur=atof(val().c_str());
        else if (a=="-ti") o.tint=atof(val().c_str());
        else if (a=="-j") o.job=val();
        else if (a=="-v") o.verbose=atoi(val().c_str());
        else if (a=="--sat-csv") o.satcsv=true;
        else if (a=="--datum") o.datum=val();
        else if (a=="-h"||a=="--help") {fputs(usage,stdout); exit(0);}
        else {fprintf(stderr,"unknown option %s\n%s",a.c_str(),usage); return false;}
    }
    if (o.conf.empty()) {fputs(usage,stderr); return false;}
    return true;
}

inline bool load_config(const CliOpts &o, nrtk_config_t &cfg, nrtk_atx_t &atx)
{
    char err[2048]="";
    nrtk_log_level(o.verbose);
    nrtk_config_default(&cfg);
    if (!nrtk_config_load(o.conf.c_str(),&cfg,err,sizeof(err))) {
        fprintf(stderr,"config error: %s\n",err); return false;
    }
    if (!o.outdir.empty()) snprintf(cfg.outdir,sizeof(cfg.outdir),"%s",o.outdir.c_str());
    if (!o.ts.empty()&&!nrtk_parse_time(o.ts.c_str(),&cfg.ts)) {fprintf(stderr,"bad -ts\n"); return false;}
    if (!o.te.empty()&&!nrtk_parse_time(o.te.c_str(),&cfg.te)) {fprintf(stderr,"bad -te\n"); return false;}
    if (o.dur>0.0&&cfg.ts.time) cfg.te=timeadd(cfg.ts,o.dur);
    if (o.tint>0.0) cfg.tint=o.tint;
    if (nrtk_config_validate(&cfg,err,sizeof(err))) {
        fprintf(stderr,"config invalid: %s\n",err); return false;
    }
    memset(&atx,0,sizeof(atx));
    if (cfg.antex[0]&&!nrtk_atx_read(cfg.antex,&atx)) {
        fprintf(stderr,"cannot read ANTEX %s\n",cfg.antex); return false;
    }
    return true;
}

} // namespace nrtk
