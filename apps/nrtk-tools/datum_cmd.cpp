/*------------------------------------------------------------------------------
 * nrtk-tools datum : output datum shift fields (datum.h)
 *
 *  check  : every point's shift, leave-one-out interpolation residuals, grid
 *           completeness; difference of the shifts between two datums
 *  at     : shift, datum coordinates and TM2 of a position
 *  theory : shift from the theory chain (station velocity + IERS Helmert)
 *           against the shift from the official coordinates
 *  tie    : station coordinates in the datums from control-point ties
 *----------------------------------------------------------------------------*/
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>
#include "nrtk/datum.h"

namespace {

const char *usage=
"usage: nrtk-tools datum check  -c conf.yaml [--datum NAME] [--time \"YYYY-MM-DD hh:mm:ss\"]\n"
"       nrtk-tools datum at     -c conf.yaml --llh lat lon h | --ecef x y z [--in NAME]\n"
"                               [--time ...] [--zone 121|119]\n"
"       nrtk-tools datum theory -c conf.yaml [--datum NAME] [--velocity-field vel.csv] [--time ...]\n"
"       nrtk-tools datum tie    -c conf.yaml --station ID --ties ties.yaml\n"
"  --in NAME : the given position is in that datum (default ITRF)\n"
"  vel.csv   : 'id,ve,vn,vu' (mm/yr, per station) or 'lat,lon,ve,vn,vu' (field, IDW)\n"
"  ties.yaml : {station: ID, ties: [{id, baseline_ecef: [dx,dy,dz] | solved_ecef: [x,y,z],\n"
"               <DATUM>: {ecef|llh}, sigma_m}]}  (baseline = control point - station, fixed solution)\n";

const char *mode_str(int m)
{
    return m==NRTK_DSHIFT_TRI?"interpolated":m==NRTK_DSHIFT_EXTRAP?"extrapolated":
           m==NRTK_DSHIFT_NEAREST?"nearest":"none";
}
void enu_at(const double *pos, const double *d, double *enu)
{
    double llh[3];
    ecef2pos(pos,llh);
    ecef2enu(llh,d,enu);
}
struct Args {
    std::string conf,datum,time,in,vel,station,ties;
    double llh[3]={0,0,0},ecef[3]={0,0,0};
    bool have_llh=false,have_ecef=false;
    double zone=0.0;
};
bool parse(int argc, char **argv, Args &a)
{
    for (int i=3;i<argc;i++) {
        std::string s=argv[i];
        auto val=[&]()->const char *{ if (i+1>=argc) {fprintf(stderr,"%s needs a value\n",s.c_str()); exit(2);} return argv[++i]; };
        if (s=="-c") a.conf=val();
        else if (s=="--datum") a.datum=val();
        else if (s=="--time") a.time=val();
        else if (s=="--in") a.in=val();
        else if (s=="--velocity-field") a.vel=val();
        else if (s=="--station") a.station=val();
        else if (s=="--ties") a.ties=val();
        else if (s=="--zone") a.zone=atof(val());
        else if (s=="--llh"&&i+3<argc) {for (int k=0;k<3;k++) a.llh[k]=atof(argv[++i]); a.have_llh=true;}
        else if (s=="--ecef"&&i+3<argc) {for (int k=0;k<3;k++) a.ecef[k]=atof(argv[++i]); a.have_ecef=true;}
        else {fputs(usage,stderr); return false;}
    }
    if (a.conf.empty()) {fputs(usage,stderr); return false;}
    return true;
}
nrtk_config_t *load(const std::string &conf)
{
    nrtk_config_t *cfg=new nrtk_config_t;
    char err[2048];
    nrtk_config_default(cfg);
    if (!nrtk_config_load(conf.c_str(),cfg,err,sizeof(err))) {fprintf(stderr,"config error: %s\n",err); return nullptr;}
    if (nrtk_config_validate(cfg,err,sizeof(err))) {fprintf(stderr,"config invalid: %s\n",err); return nullptr;}
    if (cfg->ndatum<=0) {fprintf(stderr,"no network.output_datums in %s\n",conf.c_str()); return nullptr;}
    return cfg;
}
gtime_t pick_time(const Args &a, const nrtk_config_t *cfg)
{
    gtime_t t;
    if (!a.time.empty()&&nrtk_parse_time(a.time.c_str(),&t)) return t;
    if (cfg->ts.time) return cfg->ts;
    return utc2gpst(timeget());
}
std::vector<int> datum_list(const Args &a, const nrtk_config_t *cfg)
{
    std::vector<int> v;
    if (a.datum.empty()) {for (int i=0;i<cfg->ndatum;i++) v.push_back(i); return v;}
    int k=nrtk_config_datum(cfg,a.datum.c_str());
    if (k<0) fprintf(stderr,"unknown datum %s\n",a.datum.c_str()); else v.push_back(k);
    return v;
}
void print_header(const nrtk_datum_t *d, const nrtk_config_t *cfg, gtime_t t)
{
    const nrtk_datum_def_t *def=nrtk_datum_def(d);
    char ts[40]; time2str(t,ts,0);
    int nsta=0;
    for (int i=0;i<nrtk_datum_npoints(d);i++) if (nrtk_datum_point_station(d,i)>=0) nsta++;
    printf("datum %s (%s, %s): theory %s @ %.1f; shifts at %s GPST\n",def->name,def->label,
           def->version[0]?def->version:"no version",def->frame[0]?def->frame:"?",def->epoch,ts);
    printf("  points %d (%d stations, %d calibration), triangles %d, interp %s, extrapolation limit %.0f km\n",
           nrtk_datum_npoints(d),nsta,nrtk_datum_npoints(d)-nsta,nrtk_datum_ntri(d),
           cfg->datum_interp==NRTK_DATUM_INTERP_NEAREST?"nearest":"barycentric",cfg->datum_extrap_max_km);
    if (nrtk_datum_nmissing(d)) {
        printf("  enabled stations without %s coordinates:",def->name);
        for (int k=0;k<nrtk_datum_nmissing(d);k++) printf(" %s",cfg->sta[nrtk_datum_missing(d,k)].id);
        printf("\n");
        for (int k=0;k<nrtk_datum_nmissing(d);k++) {
            printf("    %s: %s\n",cfg->sta[nrtk_datum_missing(d,k)].id,nrtk_datum_missing_covered(d,k)?
                   "covered by a calibration point":"NOT covered -> its cells are datum-incomplete (service refused there)");
        }
    }
}

/* check -----------------------------------------------------------------------------*/
int cmd_check(const Args &a)
{
    nrtk_config_t *cfg=load(a.conf);
    if (!cfg) return 1;
    gtime_t t=pick_time(a,cfg);
    std::vector<nrtk_datum_t *> ds;
    for (int i=0;i<cfg->ndatum;i++) ds.push_back(nrtk_datum_from_config(cfg,i));
    int rc=0;
    for (int k:datum_list(a,cfg)) {
        const nrtk_datum_t *d=ds[k];
        print_header(d,cfg,t);
        printf("  %-10s %-13s %8s %8s %8s | %7s %7s %7s  %s\n","id","source","dE(m)","dN(m)","dU(m)",
               "LOO dE","dN","dU(mm)","interpolation");
        double s2[3]={0,0,0}; int n=0,nbad=0;
        for (int i=0;i<nrtk_datum_npoints(d);i++) {
            double pos[3],del[3],enu[3],li[3],lenu[3]={0,0,0};
            nrtk_datum_info_t info;
            nrtk_datum_point(d,i,t,pos,del);
            enu_at(pos,del,enu);
            int mode=nrtk_datum_loo(d,i,t,li,&info);
            if (mode) {
                double r[3];
                for (int j=0;j<3;j++) r[j]=li[j]-del[j];
                enu_at(pos,r,lenu);
            }
            printf("  %-10s %-13s %8.3f %8.3f %8.3f |",nrtk_datum_point_id(d,i),nrtk_datum_point_source(d,i),
                   enu[0],enu[1],enu[2]);
            if (mode) {
                printf(" %7.1f %7.1f %7.1f  %s",lenu[0]*1E3,lenu[1]*1E3,lenu[2]*1E3,mode_str(mode));
                if (mode==NRTK_DSHIFT_TRI) {
                    for (int j=0;j<3;j++) s2[j]+=lenu[j]*lenu[j];
                    n++;
                    if (hypot(lenu[0],lenu[1])>0.05||fabs(lenu[2])>0.05) {nbad++; printf("  <-- check ITRF coordinates");}
                }
                else printf(" (hull point, not a test)");
            }
            else printf(" %7s %7s %7s  -","","","");
            printf("\n");
        }
        if (n) {
            printf("  leave-one-out residual RMS of %d interior points: E %.1f N %.1f U %.1f mm; > 50 mm: %d\n",
                   n,sqrt(s2[0]/n)*1E3,sqrt(s2[1]/n)*1E3,sqrt(s2[2]/n)*1E3,nbad);
            if (nbad) rc=3;
        }
        else printf("  no interior point for a leave-one-out test (need at least 4 points)\n");
        printf("\n");
    }
    /* difference of the shifts between datums at common points */
    for (int k1=0;k1<cfg->ndatum;k1++) for (int k2=k1+1;k2<cfg->ndatum;k2++) {
        if (!a.datum.empty()) break;
        printf("delta(%s) - delta(%s) at common points (E N U, m): should match the difference of\n"
               "the two published coordinate sets of nearby control points (Taiwan 2010/2020: 0.2-0.4 m)\n",
               cfg->datum[k2].name,cfg->datum[k1].name);
        double s[3]={0,0,0},s2[3]={0,0,0}; int n=0;
        for (int i=0;i<nrtk_datum_npoints(ds[k1]);i++) {
            const char *id=nrtk_datum_point_id(ds[k1],i);
            int j=-1;
            for (int q=0;q<nrtk_datum_npoints(ds[k2]);q++) if (!strcmp(nrtk_datum_point_id(ds[k2],q),id)) {j=q; break;}
            if (j<0) continue;
            double p[3],d1[3],d2[3],dd[3],enu[3];
            nrtk_datum_point(ds[k1],i,t,p,d1);
            nrtk_datum_point(ds[k2],j,t,nullptr,d2);
            for (int q=0;q<3;q++) dd[q]=d2[q]-d1[q];
            enu_at(p,dd,enu);
            printf("  %-10s %8.3f %8.3f %8.3f\n",id,enu[0],enu[1],enu[2]);
            for (int q=0;q<3;q++) {s[q]+=enu[q]; s2[q]+=enu[q]*enu[q];}
            n++;
        }
        if (n) {
            printf("  mean   ");
            for (int q=0;q<3;q++) printf(" %8.3f",s[q]/n);
            printf("   std  ");
            for (int q=0;q<3;q++) printf(" %8.3f",sqrt(std::max(0.0,s2[q]/n-(s[q]/n)*(s[q]/n))));
            printf("  (n=%d)\n",n);
        }
        else printf("  no common points\n");
    }
    for (auto d:ds) nrtk_datum_free(d);
    delete cfg;
    return rc;
}

/* at --------------------------------------------------------------------------------*/
int cmd_at(const Args &a)
{
    if (!a.have_llh&&!a.have_ecef) {fputs(usage,stderr); return 2;}
    nrtk_config_t *cfg=load(a.conf);
    if (!cfg) return 1;
    gtime_t t=pick_time(a,cfg);
    double pin[3],pos[3];
    if (a.have_ecef) memcpy(pin,a.ecef,sizeof(pin));
    else {double l[3]={a.llh[0]*D2R,a.llh[1]*D2R,a.llh[2]}; pos2ecef(l,pin);}
    std::vector<nrtk_datum_t *> ds;
    for (int i=0;i<cfg->ndatum;i++) ds.push_back(nrtk_datum_from_config(cfg,i));
    memcpy(pos,pin,sizeof(pos));
    if (!a.in.empty()) {       /* given in a datum: back to ITRF first */
        int k=nrtk_config_datum(cfg,a.in.c_str());
        if (k<0) {fprintf(stderr,"unknown datum %s\n",a.in.c_str()); return 2;}
        nrtk_datum_inverse(ds[k],pin,t,pos);
    }
    double llh[3],E,N;
    ecef2pos(pos,llh);
    double zone=a.zone>0.0?a.zone:nrtk_tm2_zone(llh[1]);
    nrtk_tm2(llh,zone,&E,&N);
    char ts[40]; time2str(t,ts,0);
    printf("ITRF (%s) at %s: lat %.9f lon %.9f h %.4f  ecef %.4f %.4f %.4f\n",cfg->frame,ts,
           llh[0]*R2D,llh[1]*R2D,llh[2],pos[0],pos[1],pos[2]);
    printf("  TM2 zone %.0f: E %.4f N %.4f (projection of the ITRF coordinates, for reference only)\n",zone,E,N);
    for (int k:datum_list(a,cfg)) {
        double del[3],pd[3],ld[3];
        nrtk_datum_info_t info;
        int mode=nrtk_datum_shift(ds[k],pos,t,del,&info);
        for (int q=0;q<3;q++) pd[q]=pos[q]+del[q];
        ecef2pos(pd,ld);
        nrtk_tm2(ld,zone,&E,&N);
        printf("%s (%s): delta E %.4f N %.4f U %.4f m, %s",cfg->datum[k].name,cfg->datum[k].label,
               info.enu[0],info.enu[1],info.enu[2],mode_str(mode));
        if (mode==NRTK_DSHIFT_EXTRAP||mode==NRTK_DSHIFT_NEAREST) printf(" (%.1f km outside the points)",info.extrap_km);
        if (info.theory) printf(" [theory]");
        printf("\n  points:");
        for (int q=0;q<3;q++) if (info.pt[q]>=0) printf(" %s(%.3f)",nrtk_datum_point_id(ds[k],info.pt[q]),info.w[q]);
        printf("\n  lat %.9f lon %.9f h %.4f  ecef %.4f %.4f %.4f\n  TM2 zone %.0f: E %.4f N %.4f\n",
               ld[0]*R2D,ld[1]*R2D,ld[2],pd[0],pd[1],pd[2],zone,E,N);
    }
    for (auto d:ds) nrtk_datum_free(d);
    delete cfg;
    return 0;
}

/* theory ---------------------------------------------------------------------------*/
struct VelField {
    bool loaded=false;
    std::map<std::string,std::array<double,3>> byid;              /* mm/yr ENU */
    std::vector<std::array<double,5>> grid;                       /* lat lon ve vn vu */
    bool get(const char *id, const double *pos, double *venu) const {
        auto it=byid.find(id);
        if (it!=byid.end()) {for (int k=0;k<3;k++) venu[k]=it->second[k]*1E-3; return true;}
        if (grid.empty()) return false;
        double llh[3]; ecef2pos(pos,llh);
        std::vector<std::pair<double,int>> d;
        for (size_t i=0;i<grid.size();i++) {
            double dlat=(grid[i][0]-llh[0]*R2D)*111.0,dlon=(grid[i][1]-llh[1]*R2D)*111.0*cos(llh[0]);
            d.push_back({hypot(dlat,dlon),(int)i});
        }
        std::sort(d.begin(),d.end());
        double sw=0; venu[0]=venu[1]=venu[2]=0.0;
        for (size_t i=0;i<d.size()&&i<6;i++) {
            double w=1.0/std::max(d[i].first*d[i].first,1E-6);
            for (int k=0;k<3;k++) venu[k]+=w*grid[d[i].second][2+k]*1E-3;
            sw+=w;
        }
        for (int k=0;k<3;k++) venu[k]/=sw;
        return true;
    }
};
bool read_velfield(const std::string &file, VelField &vf)
{
    std::ifstream f(file);
    if (!f) {fprintf(stderr,"cannot read %s\n",file.c_str()); return false;}
    std::string line;
    while (std::getline(f,line)) {
        for (auto &c:line) if (c==','||c==';'||c=='\t') c=' ';
        std::istringstream is(line);
        std::vector<std::string> tok;
        std::string s;
        while (is>>s) tok.push_back(s);
        if (tok.empty()||tok[0][0]=='#') continue;
        char *end;
        strtod(tok[0].c_str(),&end);
        bool numeric=*end=='\0';
        if (!numeric) {
            if (tok.size()<4) continue;         /* header line or id without values */
            strtod(tok[1].c_str(),&end);
            if (*end) continue;
            vf.byid[tok[0]]={atof(tok[1].c_str()),atof(tok[2].c_str()),atof(tok[3].c_str())};
        }
        else if (tok.size()>=5) {
            vf.grid.push_back({atof(tok[0].c_str()),atof(tok[1].c_str()),atof(tok[2].c_str()),
                               atof(tok[3].c_str()),atof(tok[4].c_str())});
        }
    }
    vf.loaded=true;
    return true;
}
int cmd_theory(const Args &a)
{
    nrtk_config_t *cfg=load(a.conf);
    if (!cfg) return 1;
    gtime_t t=pick_time(a,cfg);
    double ty=nrtk_decimal_year(t);
    VelField vf;
    if (!a.vel.empty()&&!read_velfield(a.vel,vf)) return 1;
    int rc=0;
    for (int k:datum_list(a,cfg)) {
        const nrtk_datum_def_t &def=cfg->datum[k];
        nrtk_datum_t *d=nrtk_datum_from_config(cfg,k);
        print_header(d,cfg,t);
        if (!def.frame[0]||def.epoch<=0.0) {
            printf("  no theory frame/epoch for this datum (output_datums[].frame/epoch)\n\n");
            nrtk_datum_free(d);
            continue;
        }
        printf("  method A: official - ITRF(t); method B: ITRF2020(t) moved back to %.1f with the\n"
               "  station velocity, then IERS ITRF2020 -> %s; coseismic steps are not modelled\n",def.epoch,def.frame);
        printf("  %-10s %7s | %8s %8s %8s | %8s %8s %8s | %7s %7s %7s\n","id","vel",
               "A dE","dN","dU","B dE","dN","dU","A-B dE","dN","dU(mm)");
        int nbad=0,nnovel=0;
        for (int i=0;i<nrtk_datum_npoints(d);i++) {
            const char *id=nrtk_datum_point_id(d,i);
            int s=nrtk_datum_point_station(d,i);
            double pos0[3],vel[3]={0,0,0},epoch0=0.0;
            if (s>=0) {memcpy(pos0,cfg->sta[s].pos,sizeof(pos0)); memcpy(vel,cfg->sta[s].vel,sizeof(vel)); epoch0=cfg->sta[s].coord_epoch;}
            else {
                for (int q=0;q<cfg->ndpoint;q++) if (!strcmp(cfg->dpoint[q].id,id)) {
                    memcpy(pos0,cfg->dpoint[q].pos,sizeof(pos0)); memcpy(vel,cfg->dpoint[q].vel,sizeof(vel));
                    epoch0=cfg->dpoint[q].coord_epoch;
                }
            }
            const char *vsrc="config";
            double venu[3];
            if (vf.loaded&&vf.get(id,pos0,venu)) {
                double llh[3]; ecef2pos(pos0,llh); enu2ecef(llh,venu,vel); vsrc="field";
            }
            if (norm(vel,3)<1E-6) {vsrc="none"; nnovel++;}
            double post[3],dA[3],xB[3],dB[3],eA[3],eB[3],dd[3],ed[3];
            nrtk_datum_point(d,i,t,post,dA);
            if (!nrtk_datum_theory(pos0,epoch0>0.0?epoch0:ty,vel,def.frame,def.epoch,xB)) {
                printf("  unknown frame %s\n",def.frame); break;
            }
            for (int q=0;q<3;q++) {dB[q]=xB[q]-post[q]; dd[q]=dA[q]-dB[q];}
            enu_at(post,dA,eA); enu_at(post,dB,eB); enu_at(post,dd,ed);
            bool bad=hypot(ed[0],ed[1])>0.05||fabs(ed[2])>0.08;
            if (bad&&strcmp(vsrc,"none")) nbad++;
            printf("  %-10s %7s | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f | %7.0f %7.0f %7.0f%s\n",id,vsrc,
                   eA[0],eA[1],eA[2],eB[0],eB[1],eB[2],ed[0]*1E3,ed[1]*1E3,ed[2]*1E3,
                   bad&&strcmp(vsrc,"none")?"  <-- check":"");
        }
        if (nnovel) printf("  %d points without velocity: their A-B contains %.1f years of station motion\n",
                           nnovel,ty-def.epoch);
        if (nbad) {printf("  %d points differ from the theory chain by more than 5 cm (8 cm up)\n",nbad); rc=3;}
        printf("\n");
        nrtk_datum_free(d);
    }
    delete cfg;
    return rc;
}

/* tie -------------------------------------------------------------------------------*/
int cmd_tie(const Args &a)
{
    if (a.ties.empty()) {fputs(usage,stderr); return 2;}
    nrtk_config_t *cfg=load(a.conf);
    if (!cfg) return 1;
    YAML::Node y;
    try {y=YAML::LoadFile(a.ties);} catch (const std::exception &e) {fprintf(stderr,"%s: %s\n",a.ties.c_str(),e.what()); return 1;}
    std::string staid=a.station;
    if (staid.empty()&&y["station"]) staid=y["station"].as<std::string>();
    int is=nrtk_config_station(cfg,staid.c_str());
    if (is<0) {fprintf(stderr,"unknown station %s\n",staid.c_str()); return 2;}
    const nrtk_station_t &st=cfg->sta[is];
    YAML::Node ties=y["ties"];
    if (!ties||!ties.IsSequence()) {fprintf(stderr,"%s: ties sequence required\n",a.ties.c_str()); return 1;}

    struct Tie { std::string id; double bl[3]; double sigma; std::map<int,std::array<double,3>> off; };
    std::vector<Tie> tv;
    for (auto n:ties) {
        Tie tt;
        tt.id=n["id"]?n["id"].as<std::string>():"?";
        tt.sigma=n["sigma_m"]?n["sigma_m"].as<double>():0.01;
        if (n["baseline_ecef"]) for (int k=0;k<3;k++) tt.bl[k]=n["baseline_ecef"][k].as<double>();
        else if (n["solved_ecef"]) for (int k=0;k<3;k++) tt.bl[k]=n["solved_ecef"][k].as<double>()-st.pos[k];
        else {fprintf(stderr,"tie %s: baseline_ecef or solved_ecef required\n",tt.id.c_str()); return 1;}
        for (int k=0;k<cfg->ndatum;k++) {
            YAML::Node dn=n[cfg->datum[k].name];
            if (!dn) continue;
            std::array<double,3> p{};
            if (dn["ecef"]) for (int q=0;q<3;q++) p[q]=dn["ecef"][q].as<double>();
            else if (dn["llh"]) {
                double l[3]={dn["llh"][0].as<double>()*D2R,dn["llh"][1].as<double>()*D2R,dn["llh"][2].as<double>()};
                pos2ecef(l,p.data());
            }
            else continue;
            tt.off[k]=p;
        }
        tv.push_back(tt);
    }
    printf("station %s: ITRF %.4f %.4f %.4f (epoch %.2f), %zu ties\n",st.id,st.pos[0],st.pos[1],st.pos[2],
           st.coord_epoch,tv.size());
    std::string yaml="    datum:\n";
    for (int k=0;k<cfg->ndatum;k++) {
        std::vector<int> use;
        for (size_t i=0;i<tv.size();i++) if (tv[i].off.count(k)) use.push_back((int)i);
        if (use.empty()) continue;
        printf("%s (%s): station coordinates = control point - baseline\n",cfg->datum[k].name,cfg->datum[k].label);
        /* first estimate: component-wise median (a mean of three ties is pulled
           1/3 of the way towards an outlier and would hide it), reject ties more
           than 3 cm from it, then the weighted mean of the accepted ones */
        double mean[3]={0,0,0};
        std::vector<int> ok(tv.size(),1);
        for (int q=0;q<3;q++) {
            std::vector<double> v;
            for (int i:use) v.push_back(tv[i].off[k][q]-tv[i].bl[q]);
            std::sort(v.begin(),v.end());
            mean[q]=v.size()%2?v[v.size()/2]:0.5*(v[v.size()/2-1]+v[v.size()/2]);
        }
        for (int i:use) {
            double r[3],enu[3];
            for (int q=0;q<3;q++) r[q]=(tv[i].off[k][q]-tv[i].bl[q])-mean[q];
            enu_at(st.pos,r,enu);
            if (hypot(enu[0],enu[1])>0.03||fabs(enu[2])>0.03) ok[i]=0;
        }
        {
            double sw=0,m[3]={0,0,0};
            for (int i:use) {
                if (!ok[i]) continue;
                double w=1.0/(tv[i].sigma*tv[i].sigma);
                for (int q=0;q<3;q++) m[q]+=w*(tv[i].off[k][q]-tv[i].bl[q]);
                sw+=w;
            }
            if (sw>0.0) for (int q=0;q<3;q++) mean[q]=m[q]/sw;
        }
        double s2=0; int n=0;
        for (int i:use) {
            double r[3],enu[3];
            for (int q=0;q<3;q++) r[q]=(tv[i].off[k][q]-tv[i].bl[q])-mean[q];
            enu_at(st.pos,r,enu);
            printf("  %-10s residual E %7.1f N %7.1f U %7.1f mm%s\n",tv[i].id.c_str(),enu[0]*1E3,enu[1]*1E3,enu[2]*1E3,
                   ok[i]?"":"  <-- rejected (> 3 cm)");
            if (ok[i]) {s2+=enu[0]*enu[0]+enu[1]*enu[1]+enu[2]*enu[2]; n++;}
        }
        if (!n) {printf("  no accepted tie\n"); continue;}
        double sigma=std::max(0.01,sqrt(s2/n/3.0)/sqrt((double)n));
        double delta[3],enu[3];
        for (int q=0;q<3;q++) delta[q]=mean[q]-st.pos[q];
        enu_at(st.pos,delta,enu);
        printf("  -> %.4f %.4f %.4f  (shift from ITRF E %.3f N %.3f U %.3f m, %d ties, sigma %.3f m)\n",
               mean[0],mean[1],mean[2],enu[0],enu[1],enu[2],n,sigma);
        char b[256];
        snprintf(b,sizeof(b),"      %s: {ecef: [%.4f, %.4f, %.4f], source: control-tie, sigma_m: %.3f}\n",
                 cfg->datum[k].name,mean[0],mean[1],mean[2],sigma);
        yaml+=b;
    }
    printf("\nYAML for station %s (add 'observed: <date>' and keep the tie files with the configuration):\n%s",
           st.id,yaml.c_str());
    delete cfg;
    return 0;
}

} // namespace

int datum_cmd(int argc, char **argv)
{
    if (argc<3) {fputs(usage,stderr); return 2;}
    std::string sub=argv[2];
    Args a;
    if (!parse(argc,argv,a)) return 2;
    if (sub=="check") return cmd_check(a);
    if (sub=="at") return cmd_at(a);
    if (sub=="theory") return cmd_theory(a);
    if (sub=="tie") return cmd_tie(a);
    fputs(usage,stderr);
    return 2;
}
