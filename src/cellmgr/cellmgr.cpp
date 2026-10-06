/*------------------------------------------------------------------------------
 * cellmgr.cpp : network cells - Delaunay triangulation of healthy stations,
 *               master/auxiliary selection for a virtual position
 *----------------------------------------------------------------------------*/
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <mutex>
#include <vector>
#include "nrtk/cellmgr.h"

namespace {

struct Tri { int v[3]; double cx,cy,r2; };

bool circum(const double *x, const double *y, int a, int b, int c, Tri &t)
{
    double ax=x[a],ay=y[a],bx=x[b],by=y[b],cx=x[c],cy=y[c];
    double d=2.0*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by));
    if (fabs(d)<1E-12) return false;
    double a2=ax*ax+ay*ay,b2=bx*bx+by*by,c2=cx*cx+cy*cy;
    t.cx=(a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d;
    t.cy=(a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d;
    t.r2=(ax-t.cx)*(ax-t.cx)+(ay-t.cy)*(ay-t.cy);
    /* counter-clockwise order */
    if ((bx-ax)*(cy-ay)-(by-ay)*(cx-ax)<0) {t.v[0]=a; t.v[1]=c; t.v[2]=b;}
    else {t.v[0]=a; t.v[1]=b; t.v[2]=c;}
    return true;
}

} // namespace

/* Bowyer-Watson Delaunay triangulation -------------------------------------*/
extern "C" int nrtk_delaunay(const double *xin, const double *yin, int n,
                             int (*tri)[3], int max)
{
    if (n<3) return 0;
    std::vector<double> x(xin,xin+n), y(yin,yin+n);
    double xmin=*std::min_element(x.begin(),x.end()), xmax=*std::max_element(x.begin(),x.end());
    double ymin=*std::min_element(y.begin(),y.end()), ymax=*std::max_element(y.begin(),y.end());
    double d=std::max(xmax-xmin,ymax-ymin)*20.0+1.0, mx=(xmin+xmax)/2, my=(ymin+ymax)/2;
    x.push_back(mx-d); y.push_back(my-d);
    x.push_back(mx+d); y.push_back(my-d);
    x.push_back(mx);   y.push_back(my+d);
    std::vector<Tri> tris;
    Tri t0;
    circum(x.data(),y.data(),n,n+1,n+2,t0);
    tris.push_back(t0);

    for (int p=0;p<n;p++) {
        std::vector<std::pair<int,int>> edges;
        std::vector<Tri> keep;
        for (auto &t:tris) {
            double dx=x[p]-t.cx, dy=y[p]-t.cy;
            if (dx*dx+dy*dy<t.r2*(1.0+1E-12)) {
                for (int k=0;k<3;k++) edges.push_back({t.v[k],t.v[(k+1)%3]});
            }
            else keep.push_back(t);
        }
        /* boundary of the polygonal hole: edges not shared */
        for (size_t i=0;i<edges.size();i++) {
            bool shared=false;
            for (size_t j=0;j<edges.size();j++) {
                if (i!=j&&edges[i].first==edges[j].second&&edges[i].second==edges[j].first) {
                    shared=true; break;
                }
            }
            if (shared) continue;
            Tri t;
            if (circum(x.data(),y.data(),edges[i].first,edges[i].second,p,t)) keep.push_back(t);
        }
        tris.swap(keep);
    }
    int nt=0;
    for (auto &t:tris) {
        if (t.v[0]>=n||t.v[1]>=n||t.v[2]>=n) continue;
        if (nt>=max) break;
        for (int k=0;k<3;k++) tri[nt][k]=t.v[k];
        nt++;
    }
    return nt;
}

struct nrtk_cellmgr {
    nrtk_cellmgr_opt_t opt;
    std::vector<nrtk_station_t> sta;
    std::vector<int> healthy;
    std::vector<double> e,n;            /* local plane (km) */
    double org[3];
    std::vector<std::array<int,3>> tris;
    mutable std::mutex mtx;
};

extern "C" void nrtk_cellmgr_opt_default(nrtk_cellmgr_opt_t *o)
{
    o->max_dist_km=70.0;
    o->lsm_radius_km=50.0;
    o->max_aux=6;
    o->max_baseline_km=70.0;
}

static void project(const nrtk_cellmgr_t *cm, const double *r, double *e, double *n)
{
    double dn,de;
    nrtk_ne(cm->org,r,&dn,&de);
    *e=de/1E3; *n=dn/1E3;
}
static void triangulate(nrtk_cellmgr_t *cm)
{
    std::vector<int> idx;
    std::vector<double> x,y;
    for (size_t i=0;i<cm->sta.size();i++) {
        if (!cm->healthy[i]||!cm->sta[i].enabled) continue;
        idx.push_back((int)i); x.push_back(cm->e[i]); y.push_back(cm->n[i]);
    }
    cm->tris.clear();
    if (idx.size()<3) return;
    int (*t)[3]=new int[idx.size()*3][3];
    int nt=nrtk_delaunay(x.data(),y.data(),(int)idx.size(),t,(int)idx.size()*3);
    for (int i=0;i<nt;i++) {
        /* drop degenerate slivers (all three stations nearly collinear) */
        double ax=x[t[i][1]]-x[t[i][0]],ay=y[t[i][1]]-y[t[i][0]];
        double bx=x[t[i][2]]-x[t[i][0]],by=y[t[i][2]]-y[t[i][0]];
        if (fabs(ax*by-ay*bx)<1E-3) continue;
        cm->tris.push_back({idx[t[i][0]],idx[t[i][1]],idx[t[i][2]]});
    }
    delete[] t;
}

extern "C" nrtk_cellmgr_t *nrtk_cellmgr_create(const nrtk_cellmgr_opt_t *opt,
                                               const nrtk_station_t *sta, int nsta)
{
    nrtk_cellmgr_t *cm=new nrtk_cellmgr_t;
    cm->opt=*opt;
    cm->sta.assign(sta,sta+nsta);
    cm->healthy.assign(nsta,1);
    double org[3]={0};
    int k=0;
    for (int i=0;i<nsta;i++) if (sta[i].enabled) {for (int j=0;j<3;j++) org[j]+=sta[i].pos[j]; k++;}
    for (int j=0;j<3;j++) cm->org[j]=k>0?org[j]/k:0.0;
    cm->e.resize(nsta); cm->n.resize(nsta);
    for (int i=0;i<nsta;i++) project(cm,sta[i].pos,&cm->e[i],&cm->n[i]);
    triangulate(cm);
    return cm;
}
extern "C" void nrtk_cellmgr_free(nrtk_cellmgr_t *cm) { delete cm; }

extern "C" int nrtk_cellmgr_set_health(nrtk_cellmgr_t *cm, int i, int healthy)
{
    std::lock_guard<std::mutex> lk(cm->mtx);
    if (i<0||i>=(int)cm->sta.size()) return 0;
    healthy=healthy?1:0;
    if (cm->healthy[i]==healthy) return 0;
    cm->healthy[i]=healthy;
    triangulate(cm);
    nrtk_log(NRTK_LOG_INFO,"cellmgr","station %s %s, %zu triangles",cm->sta[i].id,
             healthy?"healthy":"unhealthy",cm->tris.size());
    return 1;
}
extern "C" int nrtk_cellmgr_healthy(const nrtk_cellmgr_t *cm, int i)
{
    std::lock_guard<std::mutex> lk(cm->mtx);
    return i>=0&&i<(int)cm->sta.size()&&cm->healthy[i]&&cm->sta[i].enabled;
}
extern "C" int nrtk_cellmgr_nearest(const nrtk_cellmgr_t *cm, const double *vpos,
                                    int healthy_only, double *dist_km)
{
    std::lock_guard<std::mutex> lk(cm->mtx);
    int best=-1; double dbest=1E30;
    for (size_t i=0;i<cm->sta.size();i++) {
        if (!cm->sta[i].enabled||(healthy_only&&!cm->healthy[i])) continue;
        double d=nrtk_dist(vpos,cm->sta[i].pos);
        if (d<dbest) {dbest=d; best=(int)i;}
    }
    if (dist_km) *dist_km=best>=0?dbest/1E3:0.0;
    return best;
}
extern "C" int nrtk_cellmgr_select(const nrtk_cellmgr_t *cm, const double *vpos,
                                   nrtk_cell_t *cell)
{
    double ve,vn,dnear;
    memset(cell,0,sizeof(*cell));
    cell->tri[0]=cell->tri[1]=cell->tri[2]=-1;
    cell->nearest=nrtk_cellmgr_nearest(cm,vpos,1,&dnear);

    std::lock_guard<std::mutex> lk(cm->mtx);
    project(cm,vpos,&ve,&vn);
    if (cell->nearest<0) return 0;
    if (dnear>cm->opt.max_dist_km) {cell->reject=1; return 0;}

    std::vector<int> core;
    for (auto &t:cm->tris) {
        double x0=cm->e[t[0]],y0=cm->n[t[0]],x1=cm->e[t[1]],y1=cm->n[t[1]],x2=cm->e[t[2]],y2=cm->n[t[2]];
        double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
        double l0=((y1-y2)*(ve-x2)+(x2-x1)*(vn-y2))/den;
        double l1=((y2-y0)*(ve-x2)+(x0-x2)*(vn-y2))/den;
        double l2=1.0-l0-l1;
        if (l0>=-1E-9&&l1>=-1E-9&&l2>=-1E-9) {core={t[0],t[1],t[2]}; break;}
    }
    if (core.empty()) { /* outside: nearest three healthy stations */
        std::vector<std::pair<double,int>> d;
        for (size_t i=0;i<cm->sta.size();i++) {
            if (!cm->healthy[i]||!cm->sta[i].enabled) continue;
            d.push_back({nrtk_dist(vpos,cm->sta[i].pos),(int)i});
        }
        if (d.size()<3) return 0;
        std::sort(d.begin(),d.end());
        core={d[0].second,d[1].second,d[2].second};
        cell->extrap=1;
    }
    for (int k=0;k<3;k++) cell->tri[k]=core[k];
    /* master: nearest vertex */
    std::sort(core.begin(),core.end(),[&](int a,int b){
        return nrtk_dist(vpos,cm->sta[a].pos)<nrtk_dist(vpos,cm->sta[b].pos);});
    cell->master=core[0];
    cell->dist_km=nrtk_dist(vpos,cm->sta[cell->master].pos)/1E3;
    cell->naux=0;
    for (int k=1;k<3;k++) cell->aux[cell->naux++]=core[k];
    /* extra stations within the LSM radius of V */
    std::vector<std::pair<double,int>> extra;
    for (size_t i=0;i<cm->sta.size();i++) {
        if (!cm->healthy[i]||!cm->sta[i].enabled) continue;
        if ((int)i==core[0]||(int)i==core[1]||(int)i==core[2]) continue;
        double d=nrtk_dist(vpos,cm->sta[i].pos)/1E3;
        double dm=nrtk_dist(cm->sta[cell->master].pos,cm->sta[i].pos)/1E3;
        if (d<=cm->opt.lsm_radius_km&&dm<=cm->opt.max_baseline_km) extra.push_back({d,(int)i});
    }
    std::sort(extra.begin(),extra.end());
    for (auto &e:extra) {
        if (cell->naux>=cm->opt.max_aux) break;
        cell->aux[cell->naux++]=e.second;
    }
    cell->valid=1;
    return 1;
}
extern "C" int nrtk_cellmgr_triangles(const nrtk_cellmgr_t *cm, int (*tri)[3], int max)
{
    std::lock_guard<std::mutex> lk(cm->mtx);
    int n=0;
    for (auto &t:cm->tris) {
        if (n>=max) break;
        for (int k=0;k<3;k++) tri[n][k]=t[k];
        n++;
    }
    return n;
}
extern "C" int nrtk_cellmgr_baselines(const nrtk_cellmgr_t *cm, int (*pair)[2], int max)
{
    std::lock_guard<std::mutex> lk(cm->mtx);
    int n=0;
    for (size_t i=0;i<cm->sta.size();i++) for (size_t j=i+1;j<cm->sta.size();j++) {
        if (!cm->sta[i].enabled||!cm->sta[j].enabled) continue;
        if (nrtk_dist(cm->sta[i].pos,cm->sta[j].pos)/1E3>cm->opt.max_baseline_km) continue;
        if (n>=max) return n;
        pair[n][0]=(int)i; pair[n][1]=(int)j; n++;
    }
    return n;
}
