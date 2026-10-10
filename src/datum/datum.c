/*------------------------------------------------------------------------------
 * datum.c : output datum shift fields, theory chain, TM2 projection (datum.h)
 *----------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nrtk/datum.h"
#include "nrtk/cellmgr.h"       /* nrtk_delaunay */

typedef struct {
    char   id[16];
    int    ista;                /* station index, -1: calibration point */
    double pos[3],vel[3],epoch; /* ITRF */
    double off[3];              /* official coordinates in the datum */
    double sigma;
    char   source[16];
    int    theory;
} dpt_t;

struct nrtk_datum {
    nrtk_datum_def_t def;
    int    idatum;
    int    interp;
    double extrap_max_km;
    int    np;
    dpt_t *pt;
    double org[3];              /* origin of the local plane */
    double *e,*n;               /* point coordinates in the local plane (km) */
    int    ntri;
    int  (*tri)[3];
    int    nmissing;
    int    missing[NRTK_MAXSTA];
    int    covered[NRTK_MAXSTA];
    double cover_km;
};

/* local plane ---------------------------------------------------------------*/
static void project(const nrtk_datum_t *d, const double *r, double *e, double *n)
{
    double dn,de;
    nrtk_ne(d->org,r,&dn,&de);
    *e=de/1E3; *n=dn/1E3;
}
static void build(nrtk_datum_t *d)
{
    int i,k,nt;
    double *x,*y;
    int (*t)[3];

    d->e=d->n=NULL; d->tri=NULL; d->ntri=0;
    for (k=0;k<3;k++) d->org[k]=0.0;
    if (d->np<=0) return;
    for (i=0;i<d->np;i++) for (k=0;k<3;k++) d->org[k]+=d->pt[i].pos[k]/d->np;
    d->e=(double *)malloc(sizeof(double)*d->np);
    d->n=(double *)malloc(sizeof(double)*d->np);
    for (i=0;i<d->np;i++) project(d,d->pt[i].pos,d->e+i,d->n+i);
    if (d->np<3) return;
    x=d->e; y=d->n;
    t=(int (*)[3])malloc(sizeof(int[3])*d->np*3);
    d->tri=(int (*)[3])malloc(sizeof(int[3])*d->np*3);
    nt=nrtk_delaunay(x,y,d->np,t,d->np*3);
    for (i=0;i<nt;i++) {        /* drop slivers (three nearly collinear points) */
        double ax=x[t[i][1]]-x[t[i][0]],ay=y[t[i][1]]-y[t[i][0]];
        double bx=x[t[i][2]]-x[t[i][0]],by=y[t[i][2]]-y[t[i][0]];
        if (fabs(ax*by-ay*bx)<1E-3) continue;
        for (k=0;k<3;k++) d->tri[d->ntri][k]=t[i][k];
        d->ntri++;
    }
    free(t);
}
static void add_point(nrtk_datum_t *d, const char *id, int ista, const double *pos,
                      const double *vel, double epoch, const nrtk_stadatum_t *sd)
{
    dpt_t *p=d->pt+d->np++;
    memset(p,0,sizeof(*p));
    snprintf(p->id,sizeof(p->id),"%s",id);
    p->ista=ista;
    memcpy(p->pos,pos,sizeof(p->pos));
    memcpy(p->vel,vel,sizeof(p->vel));
    p->epoch=epoch;
    memcpy(p->off,sd->pos,sizeof(p->off));
    p->sigma=sd->sigma;
    snprintf(p->source,sizeof(p->source),"%s",sd->source);
    p->theory=!strcmp(sd->source,"theory");
}

extern nrtk_datum_t *nrtk_datum_create(const nrtk_datum_def_t *def, int idatum,
                                       const nrtk_station_t *sta, int nsta,
                                       const nrtk_dpoint_t *pts, int npts,
                                       int interp, double extrap_max_km, double cover_km)
{
    nrtk_datum_t *d=(nrtk_datum_t *)calloc(1,sizeof(nrtk_datum_t));
    int i,j,nunc=0;

    d->def=*def; d->idatum=idatum; d->interp=interp;
    d->extrap_max_km=extrap_max_km>0.0?extrap_max_km:20.0;
    d->cover_km=cover_km>0.0?cover_km:10.0;
    d->pt=(dpt_t *)calloc((size_t)(nsta+npts)+1,sizeof(dpt_t));
    for (i=0;i<nsta;i++) {
        if (sta[i].datum[idatum].valid) {
            add_point(d,sta[i].id,i,sta[i].pos,sta[i].vel,sta[i].coord_epoch,sta[i].datum+idatum);
        }
        else if (sta[i].enabled&&d->nmissing<NRTK_MAXSTA) d->missing[d->nmissing++]=i;
    }
    for (i=0;i<npts;i++) {
        if (!pts[i].datum[idatum].valid) continue;
        add_point(d,pts[i].id,-1,pts[i].pos,pts[i].vel,pts[i].coord_epoch,pts[i].datum+idatum);
    }
    /* a calibration point within cover_km stands in for a station without coordinates */
    for (i=0;i<d->nmissing;i++) {
        d->covered[i]=0;
        for (j=0;j<npts;j++) {
            if (!pts[j].datum[idatum].valid) continue;
            if (nrtk_dist(sta[d->missing[i]].pos,pts[j].pos)<=d->cover_km*1E3) {d->covered[i]=1; break;}
        }
        if (!d->covered[i]) nunc++;
    }
    build(d);
    if (d->np<3) {
        nrtk_log(NRTK_LOG_WARN,"datum","%s: only %d points with coordinates, shift field "
                 "degenerates to the nearest point",def->name,d->np);
    }
    else {
        nrtk_log(NRTK_LOG_INFO,"datum","%s (%s, %s): %d points, %d triangles, %d enabled "
                 "stations without coordinates (%d not covered by a calibration point)",
                 def->name,def->label,def->version,d->np,d->ntri,d->nmissing,nunc);
    }
    return d;
}
extern nrtk_datum_t *nrtk_datum_from_config(const nrtk_config_t *cfg, int idatum)
{
    if (idatum<0||idatum>=cfg->ndatum) return NULL;
    return nrtk_datum_create(cfg->datum+idatum,idatum,cfg->sta,cfg->nsta,cfg->dpoint,
                             cfg->ndpoint,cfg->datum_interp,cfg->datum_extrap_max_km,
                             cfg->datum_cover_km);
}
extern void nrtk_datum_free(nrtk_datum_t *d)
{
    if (!d) return;
    free(d->pt); free(d->e); free(d->n); free(d->tri);
    free(d);
}
extern const nrtk_datum_def_t *nrtk_datum_def(const nrtk_datum_t *d) {return &d->def;}
extern int nrtk_datum_npoints(const nrtk_datum_t *d) {return d->np;}
extern int nrtk_datum_ntri(const nrtk_datum_t *d) {return d->ntri;}
extern const char *nrtk_datum_point_id(const nrtk_datum_t *d, int i)
{
    return i>=0&&i<d->np?d->pt[i].id:"";
}
extern int nrtk_datum_point_station(const nrtk_datum_t *d, int i)
{
    return i>=0&&i<d->np?d->pt[i].ista:-1;
}
extern const char *nrtk_datum_point_source(const nrtk_datum_t *d, int i)
{
    return i>=0&&i<d->np?d->pt[i].source:"";
}
extern int nrtk_datum_nmissing(const nrtk_datum_t *d) {return d->nmissing;}
extern int nrtk_datum_missing(const nrtk_datum_t *d, int k)
{
    return k>=0&&k<d->nmissing?d->missing[k]:-1;
}
extern int nrtk_datum_missing_covered(const nrtk_datum_t *d, int k)
{
    return k>=0&&k<d->nmissing?d->covered[k]:0;
}
extern int nrtk_datum_nuncovered(const nrtk_datum_t *d)
{
    int k,n=0;
    for (k=0;k<d->nmissing;k++) if (!d->covered[k]) n++;
    return n;
}
extern int nrtk_datum_cell_incomplete(const nrtk_datum_t *d, const int *tri)
{
    int k,q;
    if (!tri) return 0;
    for (q=0;q<3;q++) {
        if (tri[q]<0) continue;
        for (k=0;k<d->nmissing;k++) if (d->missing[k]==tri[q]&&!d->covered[k]) return 1;
    }
    return 0;
}

/* point i: ITRF position at t and shift ---------------------------------------*/
static void point_at(const dpt_t *p, gtime_t t, double *pos, double *delta)
{
    double dy=p->epoch>0.0&&t.time?nrtk_decimal_year(t)-p->epoch:0.0;
    int k;
    for (k=0;k<3;k++) {
        double x=p->pos[k]+p->vel[k]*dy;
        if (pos) pos[k]=x;
        if (delta) delta[k]=p->off[k]-x;
    }
}
extern void nrtk_datum_point(const nrtk_datum_t *d, int i, gtime_t t, double *pos_itrf,
                             double *delta)
{
    if (i<0||i>=d->np) return;
    point_at(d->pt+i,t,pos_itrf,delta);
}

/* shift field -----------------------------------------------------------------*/
static double seg_dist(double px, double py, double ax, double ay, double bx, double by)
{
    double vx=bx-ax,vy=by-ay,wx=px-ax,wy=py-ay,l2=vx*vx+vy*vy,s;
    s=l2>0.0?(wx*vx+wy*vy)/l2:0.0;
    if (s<0.0) s=0.0; else if (s>1.0) s=1.0;
    return hypot(px-(ax+s*vx),py-(ay+s*vy));
}
static void bary(const nrtk_datum_t *d, const int *t, double ve, double vn, double *l)
{
    double x0=d->e[t[0]],y0=d->n[t[0]],x1=d->e[t[1]],y1=d->n[t[1]],x2=d->e[t[2]],y2=d->n[t[2]];
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
    l[0]=((y1-y2)*(ve-x2)+(x2-x1)*(vn-y2))/den;
    l[1]=((y2-y0)*(ve-x2)+(x0-x2)*(vn-y2))/den;
    l[2]=1.0-l[0]-l[1];
}
static int nearest_point(const nrtk_datum_t *d, double ve, double vn, double *dist)
{
    int i,best=-1;
    double db=1E30;
    for (i=0;i<d->np;i++) {
        double r=hypot(ve-d->e[i],vn-d->n[i]);
        if (r<db) {db=r; best=i;}
    }
    if (dist) *dist=db;
    return best;
}
extern int nrtk_datum_shift(const nrtk_datum_t *d, const double *pos, gtime_t t,
                            double *delta, nrtk_datum_info_t *info)
{
    nrtk_datum_info_t inf;
    double ve,vn,l[3],llh[3],dk[3];
    int i,k,best=-1;

    memset(&inf,0,sizeof(inf));
    inf.pt[0]=inf.pt[1]=inf.pt[2]=-1;
    for (k=0;k<3;k++) delta[k]=0.0;
    if (d->np<=0) {
        if (info) *info=inf;
        return NRTK_DSHIFT_NONE;
    }
    project(d,pos,&ve,&vn);
    if (d->interp==NRTK_DATUM_INTERP_NEAREST||d->ntri==0) {
        inf.mode=NRTK_DSHIFT_NEAREST;
        inf.pt[0]=nearest_point(d,ve,vn,&inf.extrap_km);
        inf.w[0]=1.0;
    }
    else {
        for (i=0;i<d->ntri;i++) {
            bary(d,d->tri[i],ve,vn,l);
            if (l[0]>=-1E-9&&l[1]>=-1E-9&&l[2]>=-1E-9) {best=i; break;}
        }
        if (best>=0) {
            inf.mode=NRTK_DSHIFT_TRI;
            for (k=0;k<3;k++) {inf.pt[k]=d->tri[best][k]; inf.w[k]=l[k]<0.0?0.0:l[k];}
        }
        else {                  /* outside: nearest triangle (planar distance to its edges) */
            double db=1E30;
            for (i=0;i<d->ntri;i++) {
                const int *tr=d->tri[i];
                double r=0.0;
                for (k=0;k<3;k++) {
                    double s=seg_dist(ve,vn,d->e[tr[k]],d->n[tr[k]],d->e[tr[(k+1)%3]],d->n[tr[(k+1)%3]]);
                    if (k==0||s<r) r=s;
                }
                if (r<db) {db=r; best=i;}
            }
            inf.extrap_km=db;
            if (db>d->extrap_max_km) {
                inf.mode=NRTK_DSHIFT_NEAREST;
                inf.pt[0]=nearest_point(d,ve,vn,NULL);
                inf.w[0]=1.0;
            }
            else {
                inf.mode=NRTK_DSHIFT_EXTRAP;
                bary(d,d->tri[best],ve,vn,l);
                for (k=0;k<3;k++) {inf.pt[k]=d->tri[best][k]; inf.w[k]=l[k];}
            }
        }
    }
    for (k=0;k<3;k++) {
        if (inf.pt[k]<0) continue;
        point_at(d->pt+inf.pt[k],t,NULL,dk);
        for (i=0;i<3;i++) delta[i]+=inf.w[k]*dk[i];
        if (d->pt[inf.pt[k]].theory) inf.theory=1;
    }
    ecef2pos(pos,llh);
    ecef2enu(llh,delta,inf.enu);
    if (info) *info=inf;
    return inf.mode;
}
extern int nrtk_datum_inverse(const nrtk_datum_t *d, const double *pos_datum, gtime_t t,
                              double *pos_itrf)
{
    double delta[3];
    int iter,k,mode=NRTK_DSHIFT_NONE;
    for (k=0;k<3;k++) pos_itrf[k]=pos_datum[k];
    for (iter=0;iter<2;iter++) {
        mode=nrtk_datum_shift(d,pos_itrf,t,delta,NULL);
        for (k=0;k<3;k++) pos_itrf[k]=pos_datum[k]-delta[k];
    }
    return mode;
}
extern int nrtk_datum_loo(const nrtk_datum_t *d, int i, gtime_t t, double *delta,
                          nrtk_datum_info_t *info)
{
    nrtk_datum_t *s;
    double pos[3];
    int j,k,mode;

    if (i<0||i>=d->np||d->np<4) return NRTK_DSHIFT_NONE;
    s=(nrtk_datum_t *)calloc(1,sizeof(nrtk_datum_t));
    s->def=d->def; s->idatum=d->idatum; s->interp=d->interp; s->extrap_max_km=d->extrap_max_km;
    s->cover_km=d->cover_km;
    s->pt=(dpt_t *)calloc((size_t)d->np,sizeof(dpt_t));
    for (j=0;j<d->np;j++) if (j!=i) s->pt[s->np++]=d->pt[j];
    build(s);
    point_at(d->pt+i,t,pos,NULL);
    mode=nrtk_datum_shift(s,pos,t,delta,info);
    if (info) {     /* point indices of the reduced set -> indices of d */
        for (k=0;k<3;k++) if (info->pt[k]>=i) info->pt[k]++;
    }
    nrtk_datum_free(s);
    return mode;
}

/* theory chain -----------------------------------------------------------------
 * IERS ITRF product centre, "Transformation parameters from ITRF2020 to past
 * ITRFs" (Transfo-ITRF2020_TRFs.txt): T (mm), D (ppb), R (mas) at 2015.0 and
 * their rates per year. X_frame = X + T + D X + R x X. */
static const struct {
    const char *frame;
    double T[3],D,R[3],Td[3],Dd,Rd[3];
} HELM[]={
    {"ITRF2014",{-1.4,-0.9,  1.4},-0.42,{0,0,0   },{0.0,-0.1, 0.2},0.00,{0,0,0   }},
    {"ITRF2008",{ 0.2, 1.0,  3.3},-0.29,{0,0,0   },{0.0,-0.1, 0.1},0.03,{0,0,0   }},
    {"ITRF2005",{ 2.7, 0.1, -1.4}, 0.65,{0,0,0   },{0.3,-0.1, 0.1},0.03,{0,0,0   }},
    {"ITRF2000",{-0.2, 0.8,-34.2}, 2.25,{0,0,0   },{0.1, 0.0,-1.7},0.11,{0,0,0   }},
    {"ITRF97",  { 6.5,-3.9,-77.9}, 3.98,{0,0,0.36},{0.1,-0.6,-3.1},0.12,{0,0,0.02}},
    {"ITRF96",  { 6.5,-3.9,-77.9}, 3.98,{0,0,0.36},{0.1,-0.6,-3.1},0.12,{0,0,0.02}},
    {"ITRF94",  { 6.5,-3.9,-77.9}, 3.98,{0,0,0.36},{0.1,-0.6,-3.1},0.12,{0,0,0.02}},
};
extern int nrtk_helmert_itrf2020(const char *frame, double epoch, const double *x,
                                 double *y)
{
    const double mas=1E-3/3600.0*D2R;
    double T[3],D,R[3],dt=epoch-2015.0;
    int i,k;

    for (i=0;i<(int)(sizeof(HELM)/sizeof(HELM[0]));i++) {
        if (!strcmp(HELM[i].frame,frame)) break;
    }
    if (i>=(int)(sizeof(HELM)/sizeof(HELM[0]))) {
        if (!strcmp(frame,"ITRF2020")) {for (k=0;k<3;k++) y[k]=x[k]; return 1;}
        return 0;
    }
    for (k=0;k<3;k++) {
        T[k]=(HELM[i].T[k]+HELM[i].Td[k]*dt)*1E-3;
        R[k]=(HELM[i].R[k]+HELM[i].Rd[k]*dt)*mas;
    }
    D=(HELM[i].D+HELM[i].Dd*dt)*1E-9;
    y[0]=x[0]+T[0]+D*x[0]-R[2]*x[1]+R[1]*x[2];
    y[1]=x[1]+T[1]+R[2]*x[0]+D*x[1]-R[0]*x[2];
    y[2]=x[2]+T[2]-R[1]*x[0]+R[0]*x[1]+D*x[2];
    return 1;
}
extern int nrtk_datum_theory(const double *pos, double t0, const double *vel,
                             const char *frame, double epoch, double *out)
{
    double x[3];
    int k;
    for (k=0;k<3;k++) x[k]=pos[k]+(vel?vel[k]:0.0)*(epoch-t0);
    return nrtk_helmert_itrf2020(frame,epoch,x,out);
}

/* TM2 (transverse Mercator, Krueger series) -------------------------------------*/
extern double nrtk_tm2_zone(double lon_rad)
{
    return lon_rad*R2D<120.0?119.0:121.0;
}
extern void nrtk_tm2(const double *pos, double lon0_deg, double *E, double *N)
{
    const double a=6378137.0,f=1.0/298.257222101,k0=0.9999,FE=250000.0;
    const double n=f/(2.0-f),n2=n*n,n3=n2*n,n4=n3*n;
    const double A=a/(1.0+n)*(1.0+n2/4.0+n4/64.0);
    const double al[4]={n/2.0-2.0*n2/3.0+5.0*n3/16.0+41.0*n4/180.0,
                        13.0*n2/48.0-3.0*n3/5.0+557.0*n4/1440.0,
                        61.0*n3/240.0-103.0*n4/140.0,
                        49561.0*n4/161280.0};
    double phi=pos[0],dl=pos[1]-lon0_deg*D2R,sp=sin(phi),c=2.0*sqrt(n)/(1.0+n);
    double t=sinh(atanh(sp)-c*atanh(c*sp));
    double xi=atan2(t,cos(dl)),eta=atanh(sin(dl)/sqrt(1.0+t*t)),xs=xi,es=eta;
    int j;
    for (j=1;j<=4;j++) {
        xs+=al[j-1]*sin(2.0*j*xi)*cosh(2.0*j*eta);
        es+=al[j-1]*cos(2.0*j*xi)*sinh(2.0*j*eta);
    }
    *E=FE+k0*A*es;
    *N=k0*A*xs;
}
