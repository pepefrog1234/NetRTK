/*------------------------------------------------------------------------------
 * nrtk-gen : offline VRS generator (milestone M3). For every offline.vrs job
 *            writes RINEX 3.04, RTCM 3 and diagnostics CSV. A job with
 *            leave_out uses that station's coordinates and excludes it.
 *----------------------------------------------------------------------------*/
#include <memory>
#include <vector>
#include "nrtk/pipeline.hpp"
#include "nrtk/datum.h"
#include "../common_cli.hpp"

using namespace nrtk;

static const char *usage=
"usage: nrtk-gen -c config.yaml [options]\n"
"  -o dir        output directory (default offline.out_dir)\n"
"  -ts/-te time  processing window (GPST \"YYYY-MM-DD hh:mm:ss\")\n"
"  -d sec        duration, -ti sec epoch interval\n"
"  -j name       run only this job\n"
"  --datum NAME  RTCM 1006 position in this output datum (network.output_datums);\n"
"                the RINEX header stays ITRF, <job>_datum.json records the shift\n"
"  -v level      log level 0-3\n";

/* <job>_datum.json: the shift applied to the RTCM 1006 position of a job */
static void write_datum_json(const std::string &file, const nrtk_datum_t *d, const double *pos,
                             const double *p1006, const nrtk_datum_info_t &info, gtime_t t)
{
    FILE *fp=fopen(file.c_str(),"w");
    if (!fp) return;
    const nrtk_datum_def_t *def=nrtk_datum_def(d);
    double llh[3],E,N;
    ecef2pos(p1006,llh);
    nrtk_tm2(llh,nrtk_tm2_zone(llh[1]),&E,&N);
    fprintf(fp,"{\"datum\":\"%s\",\"label\":\"%s\",\"version\":\"%s\",\"time\":\"%s\",\n"
               " \"itrf_ecef\":[%.4f,%.4f,%.4f],\"datum_ecef\":[%.4f,%.4f,%.4f],\n"
               " \"datum_llh\":[%.9f,%.9f,%.4f],\"tm2\":{\"zone\":%d,\"e\":%.4f,\"n\":%.4f},\n"
               " \"delta_ecef\":[%.4f,%.4f,%.4f],\"delta_enu\":[%.4f,%.4f,%.4f],\n"
               " \"mode\":\"%s\",\"extrap_km\":%.2f,\"source\":\"%s\",\"points\":[",
            def->name,def->label,def->version,time_str(t).c_str(),pos[0],pos[1],pos[2],
            p1006[0],p1006[1],p1006[2],llh[0]*R2D,llh[1]*R2D,llh[2],(int)nrtk_tm2_zone(llh[1]),E,N,
            p1006[0]-pos[0],p1006[1]-pos[1],p1006[2]-pos[2],info.enu[0],info.enu[1],info.enu[2],
            info.mode==NRTK_DSHIFT_TRI?"interpolated":info.mode==NRTK_DSHIFT_EXTRAP?"extrapolated":
            info.mode==NRTK_DSHIFT_NEAREST?"nearest":"none",info.extrap_km,
            info.theory?"theory":"official");
    int k0=0;
    for (int k=0;k<3;k++) {
        if (info.pt[k]<0) continue;
        fprintf(fp,"%s{\"id\":\"%s\",\"w\":%.4f}",k0++?",":"",nrtk_datum_point_id(d,info.pt[k]),info.w[k]);
    }
    fprintf(fp,"]}\n");
    fclose(fp);
}

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

    /* output datum of the RTCM 1006 position (--datum) */
    nrtk_datum_t *datum=nullptr;
    if (!o.datum.empty()&&o.datum!="none") {
        int id=nrtk_config_datum(cfg,o.datum.c_str());
        if (id<0) {fprintf(stderr,"--datum %s is not in network.output_datums\n",o.datum.c_str()); return 2;}
        datum=nrtk_datum_from_config(cfg,id);
    }
    OfflineReader rd(*cfg);
    if (!rd.ok()) {fprintf(stderr,"input error\n"); return 1;}
    Network net(*cfg,&atx);

    std::vector<std::unique_ptr<Job>> jobs;
    gtime_t t0=cfg->ts.time?cfg->ts:utc2gpst(timeget());
    for (int k=0;k<cfg->njob;k++) {
        const nrtk_vrsjob_t &j=cfg->job[k];
        if (!o.job.empty()&&o.job!=j.name) continue;
        auto jb=std::make_unique<Job>();
        jb->name=j.name;
        jb->exclude=-1;
        if (j.leave_out[0]) {
            jb->exclude=nrtk_config_station(cfg,j.leave_out);
            nrtk_sta_pos(cfg->sta+jb->exclude,t0,jb->pos);
        }
        else memcpy(jb->pos,j.pos,sizeof(jb->pos));
        int staid=j.staid?j.staid:3000+k;
        jb->stream=std::make_unique<VrsStream>(net,jb->pos,jb->exclude);
        double p1006[3],delta[3];
        const double *rpos=nullptr;
        if (datum) {
            nrtk_datum_info_t info;
            int mode=nrtk_datum_shift(datum,jb->pos,t0,delta,&info);
            for (int i=0;i<3;i++) p1006[i]=jb->pos[i]+delta[i];
            rpos=p1006;
            fprintf(stderr,"%s: 1006 position in %s, delta E %.3f N %.3f U %.3f m (%s%s)\n",
                    j.name,nrtk_datum_def(datum)->label,info.enu[0],info.enu[1],info.enu[2],
                    mode==NRTK_DSHIFT_TRI?"interpolated":mode==NRTK_DSHIFT_EXTRAP?"extrapolated":
                    mode==NRTK_DSHIFT_NEAREST?"nearest point":"no points",info.theory?", theory":"");
            write_datum_json(dir+"/"+jb->name+"_datum.json",datum,jb->pos,p1006,info,t0);
        }
        jb->writer=std::make_unique<VrsWriter>(*cfg,dir,jb->name,jb->pos,staid,&rd.nav(),rpos,
                                               datum?nrtk_datum_def(datum)->label:nullptr);
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
    nrtk_datum_free(datum);
    nrtk_atx_free(&atx);
    delete cfg;
    return 0;
}
