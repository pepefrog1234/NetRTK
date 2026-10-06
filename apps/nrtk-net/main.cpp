/*------------------------------------------------------------------------------
 * nrtk-net : network solution CLI - ambiguity fixing between the reference
 *            stations, fix-rate report and residual output (milestone M1)
 *----------------------------------------------------------------------------*/
#include <map>
#include <vector>
#include "nrtk/pipeline.hpp"
#include "../common_cli.hpp"

using namespace nrtk;

static const char *usage=
"usage: nrtk-net -c config.yaml [options]\n"
"  -o dir        output directory (default offline.out_dir)\n"
"  -ts time      start \"YYYY-MM-DD hh:mm:ss\" (GPST)\n"
"  -te time      end\n"
"  -d sec        duration\n"
"  -ti sec       epoch interval\n"
"  --sat-csv     write per-satellite residual CSV per baseline\n"
"  -v level      log level 0-3\n";

int main(int argc, char **argv)
{
    CliOpts o;
    nrtk_config_t *cfg=new nrtk_config_t;
    nrtk_atx_t atx;
    if (!parse_cli(argc,argv,o,usage)||!load_config(o,*cfg,atx)) return 2;
    std::string dir=cfg->outdir;
    mkdirs(dir);

    OfflineReader rd(*cfg);
    if (!rd.ok()) {fprintf(stderr,"input error\n"); return 1;}
    Network net(*cfg,&atx);

    std::map<std::pair<int,int>,FILE *> fbl,fsat;
    auto name=[&](int a,int b){return std::string(cfg->sta[a].id)+"-"+cfg->sta[b].id;};
    gtime_t t;
    std::vector<std::vector<obsd_t>> raw;
    long nep=0;
    double tproc=0.0,tmax=0.0;

    while (rd.next(t,raw)) {
        auto ep=net.process(t,raw,rd.nav_ptr());
        nep++;
        tproc+=ep->proc_ms; if (ep->proc_ms>tmax) tmax=ep->proc_ms;
        for (auto &kv:ep->bl) {
            const nrtk_blepoch_t &b=kv.second;
            FILE *&fp=fbl[kv.first];
            if (!fp) {
                fp=fopen((dir+"/bl_"+name(kv.first.first,kv.first.second)+".csv").c_str(),"w");
                fprintf(fp,"time,nsat,nfixed,nwl,stat,ratio,ztd_m,ztd_a,nslip\n");
            }
            fprintf(fp,"%s,%d,%d,%d,%d,%.2f,%.4f,%.4f,%d\n",time_str(t).c_str(),b.nsat,b.nfixed,
                    b.nwl,b.stat,b.ratio,b.ztd[0],b.ztd[1],b.nslip);
            if (!o.satcsv) continue;
            FILE *&fs=fsat[kv.first];
            if (!fs) {
                fs=fopen((dir+"/blsat_"+name(kv.first.first,kv.first.second)+".csv").c_str(),"w");
                fprintf(fs,"time,sat,el_m,az_m,fixed,wl,ewl,isref,newfix,res1,res2,res3,"
                           "code1,code2,dd_iono,dd_nondisp,dd_trop_model\n");
            }
            for (int i=0;i<b.ns;i++) {
                const nrtk_blsat_t &s=b.sat[i];
                char id[8]; satno2id(s.sat,id);
                fprintf(fs,"%s,%s,%.2f,%.2f,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.3f,%.3f,%.4f,%.4f,%.4f\n",
                        time_str(t).c_str(),id,s.azel_m[1]*R2D,s.azel_m[0]*R2D,s.fixed,s.wlfixed,
                        s.ewlfixed,s.isref,s.newfix,s.dd_res[0],s.dd_res[1],s.dd_res[2],
                        s.dd_code[0],s.dd_code[1],s.dd_iono,s.dd_nondisp,s.dd_trop_model);
            }
        }
        if (nep%600==0) {
            fprintf(stderr,"%s",time_str(t).c_str());
            for (auto &b:net.baseline_status()) {
                fprintf(stderr," %s:%.0f%%",name(b.a,b.b).c_str(),b.nepoch?100.0*b.nfixep/b.nepoch:0.0);
            }
            fprintf(stderr,"\n");
        }
    }
    for (auto &kv:fbl) fclose(kv.second);
    for (auto &kv:fsat) fclose(kv.second);

    /* summary (milestone M1 / gate G1: fix rate >= 95%, TTFF <= 120 s) */
    FILE *fs=fopen((dir+"/baselines.csv").c_str(),"w");
    fprintf(fs,"baseline,length_km,epochs,fixed_epochs,fix_rate,ttff_s,refix,slips\n");
    printf("\n%-14s %9s %8s %8s %8s %6s %6s  G1\n","baseline","len(km)","epochs","fixrate","TTFF(s)","refix","slips");
    for (auto &b:net.baseline_status()) {
        double fr=b.nepoch?100.0*b.nfixep/b.nepoch:0.0;
        bool g1=fr>=95.0&&b.ttff>=0.0&&b.ttff<=120.0;
        printf("%-14s %9.1f %8ld %7.1f%% %8.0f %6ld %6ld  %s\n",name(b.a,b.b).c_str(),b.length_km,
               b.nepoch,fr,b.ttff,b.nrefix,b.nslip,g1?"pass":"FAIL");
        fprintf(fs,"%s,%.3f,%ld,%ld,%.4f,%.1f,%ld,%ld\n",name(b.a,b.b).c_str(),b.length_km,
                b.nepoch,b.nfixep,fr/100.0,b.ttff,b.nrefix,b.nslip);
    }
    fclose(fs);
    printf("\nepochs %ld, processing %.1f ms/epoch (max %.1f ms)\n",nep,nep?tproc/nep:0.0,tmax);
    nrtk_atx_free(&atx);
    delete cfg;
    return 0;
}
