/*------------------------------------------------------------------------------
 * nrtk-tools : validation utilities
 *
 *  obscmp : observation-level comparison of a VRS RINEX file with the real
 *           observations of the left-out station at the same position
 *           (spec ch.10 "observation layer"). Per signal, double differences
 *           (VRS - station, satellite - reference) of phase and code are formed;
 *           the station observations are reduced to marker / null antenna with
 *           its ANTEX model, the integer ambiguity of each continuous arc is
 *           removed, and residual RMS are reported per signal and band.
 *----------------------------------------------------------------------------*/
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include "nrtk/config.h"
#include "nrtk/ingest.h"
#include "nrtk/antenna.h"
#include "nrtk/prep.h"

namespace {

const char *usage=
"usage: nrtk-tools obscmp -c config.yaml -s STATION -v vrs.obs [-o res.csv] [-e elmin_deg]\n"
"                         [-j jump_m]  (new arc on an unflagged DD jump, 0: off; default 0.05)\n"
"         compare VRS observations with the real observations of STATION\n"
"         (rinex file and antenna taken from the configuration)\n"
"       nrtk-tools datum check|at|theory|tie ...   output datum shift fields (datum.h)\n"
"       nrtk-tools confjson config.yaml\n";

struct Sample { gtime_t t; double el,ph,cd; };
struct Arc { int sat; uint8_t code; double lam; std::vector<Sample> s; };

double band_of(int sys, uint8_t code)
{
    double f=nrtk_code2freq(sys,code);
    return f;
}
std::string band_name(double f)
{
    double M=f/1E6;
    if (M>1550) return "L1/E1/B1";
    if (M>1200) return "L2/E5b/B3";
    return "L5/E5a";
}

int obscmp(int argc, char **argv)
{
    std::string conf,staid,vrsfile,outfile;
    double elmin=15.0,jump=0.05;
    for (int i=2;i<argc;i++) {
        std::string a=argv[i];
        if (a=="-c"&&i+1<argc) conf=argv[++i];
        else if (a=="-s"&&i+1<argc) staid=argv[++i];
        else if (a=="-v"&&i+1<argc) vrsfile=argv[++i];
        else if (a=="-o"&&i+1<argc) outfile=argv[++i];
        else if (a=="-e"&&i+1<argc) elmin=atof(argv[++i]);
        else if (a=="-j"&&i+1<argc) jump=atof(argv[++i]);
        else {fputs(usage,stderr); return 2;}
    }
    if (conf.empty()||staid.empty()||vrsfile.empty()) {fputs(usage,stderr); return 2;}

    nrtk_config_t *cfg=new nrtk_config_t;
    char err[1024];
    nrtk_config_default(cfg);
    if (!nrtk_config_load(conf.c_str(),cfg,err,sizeof(err))) {fprintf(stderr,"%s\n",err); return 1;}
    int ista=nrtk_config_station(cfg,staid.c_str());
    if (ista<0) {fprintf(stderr,"unknown station %s\n",staid.c_str()); return 1;}
    const nrtk_station_t *sta=cfg->sta+ista;
    nrtk_atx_t atx={0};
    const nrtk_ant_t *ant=nullptr;
    if (cfg->antex[0]&&nrtk_atx_read(cfg->antex,&atx)) ant=nrtk_atx_find(&atx,sta->antenna);
    if (!ant) fprintf(stderr,"warning: antenna %s not found, PCO/PCV not removed\n",sta->antenna);

    nav_t nav={0};
    if (!nrtk_nav_read(cfg->nav,&nav)) {fprintf(stderr,"no navigation data\n"); return 1;}
    nrtk_rnxsrc_t A,B;
    if (!nrtk_rnxsrc_open(&A,vrsfile.c_str(),1)||!nrtk_rnxsrc_open(&B,sta->obsfile,2)) return 1;

    std::vector<obsd_t> oa(MAXOBS),ob(MAXOBS);
    std::map<std::pair<int,uint8_t>,int> ref;          /* (sys,code) -> ref sat */
    std::map<std::tuple<int,uint8_t,int>,int> cur;     /* (sat,code,ref) -> arc index */
    std::vector<Arc> arcs;
    std::map<std::pair<int,uint8_t>,gtime_t> last;
    long nep=0;
    gtime_t tprev={0};
    double dt_nom=0.0;          /* epoch interval of the files (arc gap threshold scales with it) */

    for (;;) {
        gtime_t t=nrtk_rnxsrc_peek(&A);
        if (!t.time) break;
        int na=nrtk_rnxsrc_read(&A,t,oa.data(),MAXOBS);
        int nb=nrtk_rnxsrc_read(&B,t,ob.data(),MAXOBS);
        if (na<=0) continue;
        if (nb<=0) continue;
        nep++;
        if (tprev.time) {
            double d=timediff(t,tprev);
            if (d>0.0&&(dt_nom<=0.0||d<0.9*dt_nom)) dt_nom=d;
        }
        tprev=t;
        const double gap=std::max(5.0,2.5*dt_nom);
        double rr[3],pos[3];
        nrtk_sta_pos(sta,t,rr);
        nrtk_tide_pos(t,rr,rr);
        ecef2pos(rr,pos);
        /* single differences per (sat,code) */
        struct SD { double ph,cd,el,lam; int slip; };
        std::map<std::pair<int,uint8_t>,SD> sd;
        for (int i=0;i<na;i++) {
            const obsd_t *a=&oa[i],*b=nullptr;
            for (int j=0;j<nb;j++) if (ob[j].sat==a->sat) {b=&ob[j]; break;}
            if (!b) continue;
            int sys=satsys(a->sat,nullptr);
            if (nrtk_sysidx(sys)<0) continue;
            double rs[6],dts[2],rate,e[3],azel[2];
            int svh;
            double r=nrtk_satrange(t,rr,a->sat,&nav,rs,dts,&rate,&svh);
            if (r<=0.0) continue;
            for (int k=0;k<3;k++) e[k]=(rs[k]-rr[k])/r;
            satazel(pos,e,azel);
            if (azel[1]<10.0*D2R) continue;
            for (int ka=0;ka<NFREQ+NEXOBS;ka++) {
                uint8_t c=a->code[ka];
                if (c==CODE_NONE||a->L[ka]==0.0||a->P[ka]==0.0) continue;
                int kb=-1;
                for (int j=0;j<NFREQ+NEXOBS;j++) if (b->code[j]==c) {kb=j; break;}
                if (kb<0||b->L[kb]==0.0||b->P[kb]==0.0) continue;
                double f=nrtk_code2freq(sys,c),lam=CLIGHT/f;
                double dant=nrtk_ant_dr(ant,sys,f,azel,sta->antdel,1);
                SD v;
                v.ph=a->L[ka]*lam-(b->L[kb]*lam-dant);
                v.cd=a->P[ka]-(b->P[kb]-dant);
                v.el=azel[1]; v.lam=lam;
                v.slip=((a->LLI[ka]|b->LLI[kb])&LLI_SLIP)?1:0;
                sd[{a->sat,c}]=v;
            }
        }
        /* reference per (system,code): keep while available and continuous */
        std::map<std::pair<int,uint8_t>,int> best;
        for (auto &kv:sd) {
            int sys=satsys(kv.first.first,nullptr);
            auto key=std::make_pair(sys,kv.first.second);
            auto it=best.find(key);
            if (it==best.end()||kv.second.el>sd[{it->second,kv.first.second}].el) best[key]=kv.first.first;
        }
        for (auto &kv:best) {
            auto it=ref.find(kv.first);
            bool keep=it!=ref.end()&&sd.count({it->second,kv.first.second})&&
                      !sd[{it->second,kv.first.second}].slip;
            if (!keep) ref[kv.first]=kv.second;
        }
        for (auto &kv:sd) {
            int sat=kv.first.first; uint8_t c=kv.first.second;
            int sys=satsys(sat,nullptr);
            int r=ref[{sys,c}];
            if (r==sat||!sd.count({r,c})) continue;
            const SD &s=kv.second,&q=sd[{r,c}];
            auto key=std::make_tuple(sat,c,r);
            auto lk=std::make_pair(sat,c);
            bool newarc=!cur.count(key)||s.slip||q.slip||
                        (last.count(lk)&&timediff(t,last[lk])>gap);
            /* slips not flagged in either file: DD phase jump between epochs */
            if (!newarc) {
                const Arc &ar=arcs[cur[key]];
                if (jump>0.0&&!ar.s.empty()&&fabs((s.ph-q.ph)-ar.s.back().ph)>jump) newarc=true;
            }
            if (newarc) {
                Arc a; a.sat=sat; a.code=c; a.lam=s.lam;
                arcs.push_back(a);
                cur[key]=(int)arcs.size()-1;
            }
            last[lk]=t;
            arcs[cur[key]].s.push_back({t,s.el,s.ph-q.ph,s.cd-q.cd});
        }
    }
    /* remove the integer of each arc, statistics */
    struct Acc { double ph2=0,cd2=0; long n=0; double phmax=0; };
    std::map<std::string,Acc> bycode,byband;
    FILE *fp=outfile.empty()?nullptr:fopen(outfile.c_str(),"w");
    if (fp) fprintf(fp,"time,sat,code,el,phase_res,code_res\n");
    long nshort=0;
    for (auto &a:arcs) {
        if (a.s.size()<10) {nshort+=a.s.size(); continue;}
        std::vector<double> v;
        for (auto &s:a.s) v.push_back(s.ph);
        std::nth_element(v.begin(),v.begin()+v.size()/2,v.end());
        double N=floor(v[v.size()/2]/a.lam+0.5);
        int sys=satsys(a.sat,nullptr);
        const char sysc=sys==SYS_GAL?'E':(sys==SYS_CMP?'C':(sys==SYS_QZS?'J':'G'));
        std::string cname=std::string(1,sysc)+code2obs(a.code);
        std::string band=band_name(band_of(sys,a.code));
        char id[8]; satno2id(a.sat,id);
        for (auto &s:a.s) {
            double ph=s.ph-N*a.lam;
            if (fp) {
                char ts[40]; time2str(s.t,ts,1);
                fprintf(fp,"%s,%s,%s,%.2f,%.4f,%.3f\n",ts,id,code2obs(a.code),s.el*R2D,ph,s.cd);
            }
            if (s.el<elmin*D2R) continue;
            for (Acc *acc: {&bycode[cname],&byband[band]}) {
                acc->ph2+=ph*ph; acc->cd2+=s.cd*s.cd; acc->n++;
                acc->phmax=std::max(acc->phmax,fabs(ph));
            }
        }
    }
    if (fp) fclose(fp);
    printf("VRS %s vs %s: %ld epochs, elevation >= %.0f deg\n",vrsfile.c_str(),staid.c_str(),nep,elmin);
    printf("%-12s %9s %12s %12s %12s\n","signal","n","phase_rms_m","phase_max_m","code_rms_m");
    for (auto &kv:bycode) {
        const Acc &a=kv.second;
        printf("%-12s %9ld %12.4f %12.4f %12.3f\n",kv.first.c_str(),a.n,sqrt(a.ph2/a.n),a.phmax,sqrt(a.cd2/a.n));
    }
    printf("%-12s\n","band:");
    for (auto &kv:byband) {
        const Acc &a=kv.second;
        printf("%-12s %9ld %12.4f %12.4f %12.3f\n",kv.first.c_str(),a.n,sqrt(a.ph2/a.n),a.phmax,sqrt(a.cd2/a.n));
    }
    if (nshort) printf("(%ld samples in arcs shorter than 10 epochs ignored)\n",nshort);
    nrtk_rnxsrc_close(&A); nrtk_rnxsrc_close(&B);
    nrtk_nav_free(&nav); nrtk_atx_free(&atx);
    delete cfg;
    return 0;
}

} // namespace

/* confjson: parsed configuration as JSON (for the Python validation scripts) */
int confjson(int argc, char **argv)
{
    if (argc<3) {fputs("usage: nrtk-tools confjson config.yaml\n",stderr); return 2;}
    nrtk_config_t *c=new nrtk_config_t;
    char err[1024];
    nrtk_config_default(c);
    if (!nrtk_config_load(argv[2],c,err,sizeof(err))) {fprintf(stderr,"%s\n",err); return 1;}
    auto esc=[](const char *s){std::string o; for (;*s;s++) {if (*s=='"'||*s=='\\') o+='\\'; o+=*s;} return o;};
    double ep[6];
    printf("{\"name\":\"%s\",\"antex\":\"%s\",\"nav\":\"%s\",\"outdir\":\"%s\",",
           esc(c->name).c_str(),esc(c->antex).c_str(),esc(c->nav).c_str(),esc(c->outdir).c_str());
    time2epoch(c->ts,ep);
    printf("\"start\":\"%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%02.0f\",",ep[0],ep[1],ep[2],ep[3],ep[4],ep[5]);
    time2epoch(c->te,ep);
    printf("\"end\":\"%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%02.0f\",\"interval\":%.3f,",ep[0],ep[1],ep[2],ep[3],ep[4],ep[5],c->tint);
    printf("\"systems\":%d,\"stations\":[",c->netpos.sysmask);
    for (int i=0;i<c->nsta;i++) {
        const nrtk_station_t *s=c->sta+i;
        double pos[3]; ecef2pos(s->pos,pos);
        printf("%s{\"id\":\"%s\",\"ecef\":[%.4f,%.4f,%.4f],\"llh\":[%.9f,%.9f,%.4f],"
               "\"antenna\":\"%s\",\"antdel\":[%.4f,%.4f,%.4f],\"rinex\":\"%s\",\"enabled\":%d}",
               i?",":"",s->id,s->pos[0],s->pos[1],s->pos[2],pos[0]*R2D,pos[1]*R2D,pos[2],
               esc(s->antenna).c_str(),s->antdel[0],s->antdel[1],s->antdel[2],esc(s->obsfile).c_str(),s->enabled);
    }
    printf("],\"jobs\":[");
    for (int i=0;i<c->njob;i++) {
        const nrtk_vrsjob_t *j=c->job+i;
        printf("%s{\"name\":\"%s\",\"leave_out\":\"%s\",\"ecef\":[%.4f,%.4f,%.4f]}",i?",":"",
               j->name,j->leave_out,j->pos[0],j->pos[1],j->pos[2]);
    }
    printf("]}\n");
    delete c;
    return 0;
}

int datum_cmd(int argc, char **argv);     /* datum_cmd.cpp */

int main(int argc, char **argv)
{
    if (argc<2) {fputs(usage,stderr); return 2;}
    nrtk_log_level(NRTK_LOG_WARN);
    std::string cmd=argv[1];
    if (cmd=="obscmp") return obscmp(argc,argv);
    if (cmd=="confjson") return confjson(argc,argv);
    if (cmd=="datum") return datum_cmd(argc,argv);
    fputs(usage,stderr);
    return 2;
}
