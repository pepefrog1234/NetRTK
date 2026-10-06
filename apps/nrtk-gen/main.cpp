/*------------------------------------------------------------------------------
 * nrtk-gen : offline VRS generator (milestone M3). For every offline.vrs job
 *            writes RINEX 3.04, RTCM 3 and diagnostics CSV. A job with
 *            leave_out uses that station's coordinates and excludes it.
 *----------------------------------------------------------------------------*/
#include <memory>
#include <vector>
#include "nrtk/pipeline.hpp"
#include "../common_cli.hpp"

using namespace nrtk;

static const char *usage=
"usage: nrtk-gen -c config.yaml [options]\n"
"  -o dir        output directory (default offline.out_dir)\n"
"  -ts/-te time  processing window (GPST \"YYYY-MM-DD hh:mm:ss\")\n"
"  -d sec        duration, -ti sec epoch interval\n"
"  -j name       run only this job\n"
"  -v level      log level 0-3\n";

struct Job {
    std::string name;
    double pos[3];
    int exclude;
    std::unique_ptr<VrsStream> stream;
    std::unique_ptr<VrsWriter> writer;
    long nep=0,nvrs=0,nnear=0,nsat=0,ndeg=0,nreset=0;
};

int main(int argc, char **argv)
{
    CliOpts o;
    nrtk_config_t *cfg=new nrtk_config_t;
    nrtk_atx_t atx;
    if (!parse_cli(argc,argv,o,usage)||!load_config(o,*cfg,atx)) return 2;
    if (cfg->njob<=0) {fprintf(stderr,"no offline.vrs jobs in config\n"); return 2;}
    std::string dir=cfg->outdir;
    mkdirs(dir);

    OfflineReader rd(*cfg);
    if (!rd.ok()) {fprintf(stderr,"input error\n"); return 1;}
    Network net(*cfg,&atx);

    std::vector<std::unique_ptr<Job>> jobs;
    for (int k=0;k<cfg->njob;k++) {
        const nrtk_vrsjob_t &j=cfg->job[k];
        if (!o.job.empty()&&o.job!=j.name) continue;
        auto jb=std::make_unique<Job>();
        jb->name=j.name;
        jb->exclude=-1;
        if (j.leave_out[0]) {
            jb->exclude=nrtk_config_station(cfg,j.leave_out);
            gtime_t t0=cfg->ts.time?cfg->ts:utc2gpst(timeget());
            nrtk_sta_pos(cfg->sta+jb->exclude,t0,jb->pos);
        }
        else memcpy(jb->pos,j.pos,sizeof(jb->pos));
        int staid=j.staid?j.staid:3000+k;
        jb->stream=std::make_unique<VrsStream>(net,jb->pos,jb->exclude);
        jb->writer=std::make_unique<VrsWriter>(*cfg,dir,jb->name,jb->pos,staid,&rd.nav());
        jobs.push_back(std::move(jb));
    }
    gtime_t t;
    std::vector<std::vector<obsd_t>> raw;
    VrsOutput *out=new VrsOutput;
    long nep=0;
    while (rd.next(t,raw)) {
        auto ep=net.process(t,raw,rd.nav_ptr());
        for (auto &j:jobs) {
            j->stream->generate(*ep,*out);
            j->writer->write(*ep,*out);
            j->nep++;
            if (out->mode==1) {j->nvrs++; j->nsat+=out->n; j->nreset+=out->info.nlockreset;}
            if (out->mode==2) j->nnear++;
            if (out->degraded) j->ndeg++;
        }
        if (++nep%600==0) {
            fprintf(stderr,"%s",time_str(t).c_str());
            for (auto &j:jobs) fprintf(stderr," %s:%ldvrs",j->name.c_str(),j->nvrs);
            fprintf(stderr,"\n");
        }
    }
    /* network layer summary (same run) */
    {
        FILE *fs=fopen((dir+"/baselines.csv").c_str(),"w");
        fprintf(fs,"baseline,length_km,epochs,fixed_epochs,fix_rate,ttff_s,refix,slips\n");
        for (auto &b:net.baseline_status()) {
            fprintf(fs,"%s-%s,%.3f,%ld,%ld,%.4f,%.1f,%ld,%ld\n",cfg->sta[b.a].id,cfg->sta[b.b].id,
                    b.length_km,b.nepoch,b.nfixep,b.nepoch?(double)b.nfixep/b.nepoch:0.0,b.ttff,
                    b.nrefix,b.nslip);
        }
        fclose(fs);
    }
    printf("\n%-16s %8s %8s %8s %8s %9s %9s\n","job","epochs","vrs","nearest","degraded","sats/ep","resets/ep");
    for (auto &j:jobs) {
        j->writer->close();
        printf("%-16s %8ld %8ld %8ld %8ld %9.1f %9.2f\n",j->name.c_str(),j->nep,j->nvrs,j->nnear,
               j->ndeg,j->nvrs?(double)j->nsat/j->nvrs:0.0,j->nvrs?(double)j->nreset/j->nvrs:0.0);
    }
    delete out;
    nrtk_atx_free(&atx);
    delete cfg;
    return 0;
}
