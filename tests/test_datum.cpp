/* unit tests: output datum shift fields, theory chain, TM2 (datum.h) */
#include <gtest/gtest.h>
#include <cmath>
#include <cstring>
#include "nrtk/datum.h"

namespace {

nrtk_station_t sta_at(const char *id, double lat, double lon, double h)
{
    nrtk_station_t s;
    memset(&s,0,sizeof(s));
    snprintf(s.id,sizeof(s.id),"%s",id);
    double p[3]={lat*D2R,lon*D2R,h};
    pos2ecef(p,s.pos);
    s.enabled=1;
    s.coord_epoch=2026.0;
    return s;
}
/* official coordinates = ITRF + delta (local ENU, m) */
void set_official(nrtk_station_t &s, int k, double de, double dn, double du, const char *src="control-tie")
{
    double llh[3],enu[3]={de,dn,du},d[3];
    ecef2pos(s.pos,llh);
    enu2ecef(llh,enu,d);
    for (int i=0;i<3;i++) s.datum[k].pos[i]=s.pos[i]+d[i];
    s.datum[k].valid=1;
    snprintf(s.datum[k].source,sizeof(s.datum[k].source),"%s",src);
}
void enu_of(const double *pos, const double *d, double *enu)
{
    double llh[3];
    ecef2pos(pos,llh);
    ecef2enu(llh,d,enu);
}
/* a shift field linear in ECEF: delta = A + Bx*E_km + By*N_km with the ECEF
 * vectors A, Bx, By given in the ENU frame of the origin (lat0,lon0) */
struct LinField {
    double lat0,lon0;
    double a[3],bx[3],by[3];        /* ENU at the origin */
    void at(const double *pos, double *delta) const {
        double o[3]={lat0*D2R,lon0*D2R,0},org[3],dn,de,A[3],Bx[3],By[3];
        pos2ecef(o,org);
        enu2ecef(o,a,A); enu2ecef(o,bx,Bx); enu2ecef(o,by,By);
        nrtk_ne(org,pos,&dn,&de);
        for (int k=0;k<3;k++) delta[k]=A[k]+Bx[k]*de/1E3+By[k]*dn/1E3;
    }
    void apply(nrtk_station_t &s, int k) const {
        double d[3]; at(s.pos,d);
        for (int i=0;i<3;i++) s.datum[k].pos[i]=s.pos[i]+d[i];
        s.datum[k].valid=1;
        snprintf(s.datum[k].source,sizeof(s.datum[k].source),"control-tie");
    }
};
nrtk_datum_def_t def(const char *name)
{
    nrtk_datum_def_t d;
    memset(&d,0,sizeof(d));
    snprintf(d.name,sizeof(d.name),"%s",name);
    snprintf(d.label,sizeof(d.label),"%s",name);
    return d;
}
gtime_t t2026()
{
    double ep[]={2026,10,11,0,0,0};
    return epoch2time(ep);
}

/* numerical meridian arc (GRS80) for the TM2 test */
double meridian_arc(double phi)
{
    const double a=6378137.0,f=1.0/298.257222101,e2=f*(2-f);
    const int n=20000;
    double h=phi/n,s=0.0;
    for (int i=0;i<=n;i++) {
        double x=i*h,v=a*(1-e2)/pow(1-e2*sin(x)*sin(x),1.5);
        s+=v*(i==0||i==n?1:(i%2?4:2));
    }
    return s*h/3.0;
}
/* Snyder (USGS PP 1395) transverse Mercator series, independent of the Krueger form */
void snyder_tm(double phi, double lam, double lam0, double *E, double *N)
{
    const double a=6378137.0,f=1.0/298.257222101,e2=f*(2-f),ep2=e2/(1-e2),k0=0.9999;
    double sp=sin(phi),cp=cos(phi),tp=tan(phi);
    double Nn=a/sqrt(1-e2*sp*sp),T=tp*tp,C=ep2*cp*cp,A=(lam-lam0)*cp;
    double e4=e2*e2,e6=e4*e2;
    double M=a*((1-e2/4-3*e4/64-5*e6/256)*phi-(3*e2/8+3*e4/32+45*e6/1024)*sin(2*phi)
                +(15*e4/256+45*e6/1024)*sin(4*phi)-(35*e6/3072)*sin(6*phi));
    *E=250000.0+k0*Nn*(A+(1-T+C)*pow(A,3)/6+(5-18*T+T*T+72*C-58*ep2)*pow(A,5)/120);
    *N=k0*(M+Nn*tp*(A*A/2+(5-T+9*C+4*C*C)*pow(A,4)/24+(61-58*T+T*T+600*C-330*ep2)*pow(A,6)/720));
}

} // namespace

TEST(Datum, HelmertItrf2020ToItrf94At2010)
{
    /* Taiwan: the IERS parameters give north -5.5 cm, east +0.3 cm, up -0.7 cm */
    double llh[3]={24.0*D2R,121.0*D2R,100.0},x[3],y[3],d[3],enu[3];
    pos2ecef(llh,x);
    ASSERT_TRUE(nrtk_helmert_itrf2020("ITRF94",2010.0,x,y));
    for (int k=0;k<3;k++) d[k]=y[k]-x[k];
    ecef2enu(llh,d,enu);
    EXPECT_NEAR(enu[1],-0.055,0.003);
    EXPECT_NEAR(enu[0], 0.003,0.003);
    EXPECT_NEAR(enu[2],-0.007,0.004);
    /* ITRF2014 at 2020.0: millimetres */
    ASSERT_TRUE(nrtk_helmert_itrf2020("ITRF2014",2020.0,x,y));
    for (int k=0;k<3;k++) d[k]=y[k]-x[k];
    ecef2enu(llh,d,enu);
    for (int k=0;k<3;k++) EXPECT_LT(fabs(enu[k]),0.005);
    EXPECT_TRUE(nrtk_helmert_itrf2020("ITRF2020",2015.0,x,y));
    for (int k=0;k<3;k++) EXPECT_DOUBLE_EQ(y[k],x[k]);
    EXPECT_FALSE(nrtk_helmert_itrf2020("TWD97",2010.0,x,y));
    /* theory chain: back by the velocity, identity frame */
    double vel[3]={0.02,-0.01,0.005},out[3];
    ASSERT_TRUE(nrtk_datum_theory(x,2026.74,vel,"ITRF2020",2010.0,out));
    for (int k=0;k<3;k++) EXPECT_NEAR(out[k],x[k]-vel[k]*16.74,1E-9);
}

TEST(Datum, Tm2MatchesMeridianArcAndSnyder)
{
    double pos[3]={24.5*D2R,121.0*D2R,0},E,N;
    nrtk_tm2(pos,121.0,&E,&N);
    EXPECT_NEAR(E,250000.0,1E-6);
    EXPECT_NEAR(N,0.9999*meridian_arc(24.5*D2R),5E-4);
    double pts[][2]={{24.5,121.8},{25.3,121.55},{22.0,120.3},{23.6,120.2}};
    for (auto &p:pts) {
        double q[3]={p[0]*D2R,p[1]*D2R,0},Es,Ns;
        nrtk_tm2(q,121.0,&E,&N);
        snyder_tm(q[0],q[1],121.0*D2R,&Es,&Ns);
        EXPECT_NEAR(E,Es,1E-3)<<p[0]<<" "<<p[1];
        EXPECT_NEAR(N,Ns,1E-3)<<p[0]<<" "<<p[1];
    }
    /* Penghu: zone 119 */
    double ph[3]={23.57*D2R,119.58*D2R,0},Es,Ns;
    EXPECT_DOUBLE_EQ(nrtk_tm2_zone(ph[1]),119.0);
    EXPECT_DOUBLE_EQ(nrtk_tm2_zone(120.5*D2R),121.0);
    nrtk_tm2(ph,119.0,&E,&N);
    snyder_tm(ph[0],ph[1],119.0*D2R,&Es,&Ns);
    EXPECT_NEAR(E,Es,1E-3);
    EXPECT_NEAR(N,Ns,1E-3);
}

TEST(Datum, BarycentricExactForLinearFieldAndContinuous)
{
    nrtk_station_t sta[4]={sta_at("A",24.95,121.25,100),sta_at("B",24.70,120.92,50),
                           sta_at("C",24.62,121.30,600),sta_at("D",25.05,120.85,30)};
    LinField lf={24.85,121.05,{0.45,-0.30,-0.02},{0.004,0.002,0.0005},{-0.003,0.005,-0.0008}};
    for (auto &s:sta) lf.apply(s,0);
    nrtk_datum_def_t dd=def("TWD97_2010");
    nrtk_datum_t *d=nrtk_datum_create(&dd,0,sta,4,nullptr,0,NRTK_DATUM_INTERP_BARY,20.0,10.0);
    EXPECT_EQ(nrtk_datum_npoints(d),4);
    EXPECT_EQ(nrtk_datum_ntri(d),2);
    EXPECT_EQ(nrtk_datum_nmissing(d),0);
    gtime_t t=t2026();
    /* interior points: a linear field is reproduced exactly */
    double pts[][2]={{24.85,121.15},{24.80,121.05},{24.90,120.95},{24.75,121.20},{24.95,121.0}};
    for (auto &p:pts) {
        double llh[3]={p[0]*D2R,p[1]*D2R,200},pos[3],delta[3],enu[3],want[3];
        nrtk_datum_info_t info;
        pos2ecef(llh,pos);
        int mode=nrtk_datum_shift(d,pos,t,delta,&info);
        EXPECT_EQ(mode,NRTK_DSHIFT_TRI)<<p[0]<<" "<<p[1];
        EXPECT_NEAR(info.w[0]+info.w[1]+info.w[2],1.0,1E-9);
        lf.at(pos,want);
        for (int k=0;k<3;k++) EXPECT_NEAR(delta[k],want[k],2E-5)<<p[0]<<" "<<p[1]<<" k="<<k;
        enu_of(pos,delta,enu);
        for (int k=0;k<3;k++) EXPECT_NEAR(info.enu[k],enu[k],1E-9);
    }
    /* a vertex: exactly its own shift */
    double delta[3],dv[3];
    nrtk_datum_info_t info;
    EXPECT_EQ(nrtk_datum_shift(d,sta[2].pos,t,delta,&info),NRTK_DSHIFT_TRI);
    for (int k=0;k<3;k++) EXPECT_NEAR(delta[k],sta[2].datum[0].pos[k]-sta[2].pos[k],1E-6);
    nrtk_datum_point(d,2,t,nullptr,dv);
    for (int k=0;k<3;k++) EXPECT_NEAR(dv[k],delta[k],1E-6);
    nrtk_datum_free(d);

    /* non-linear field (independent random shifts): continuous across the
       shared edge of the two triangles */
    double rnd[4][3]={{0.41,-0.33,0.02},{0.52,-0.21,-0.04},{0.38,-0.40,0.07},{0.60,-0.25,0.01}};
    for (int i=0;i<4;i++) set_official(sta[i],0,rnd[i][0],rnd[i][1],rnd[i][2]);
    d=nrtk_datum_create(&dd,0,sta,4,nullptr,0,NRTK_DATUM_INTERP_BARY,20.0,10.0);
    /* find the shared edge: the two vertices present in both triangles */
    double mid[3]={0,0,0},dir[3]={0,0,0};
    {
        int tri[8][3]; int nt=0;
        /* reuse the delaunay of the same four points in the local plane */
        double org[3]={24.85*D2R,121.05*D2R,0},o[3],x[4],y[4];
        pos2ecef(org,o);
        for (int i=0;i<4;i++) {double dn,de; nrtk_ne(o,sta[i].pos,&dn,&de); x[i]=de; y[i]=dn;}
        nt=nrtk_delaunay(x,y,4,tri,8);
        ASSERT_EQ(nt,2);
        int shared[2],ns=0;
        for (int k=0;k<3;k++) for (int q=0;q<3;q++) if (tri[0][k]==tri[1][q]&&ns<2) shared[ns++]=tri[0][k];
        ASSERT_EQ(ns,2);
        for (int k=0;k<3;k++) mid[k]=0.5*(sta[shared[0]].pos[k]+sta[shared[1]].pos[k]);
        /* direction across the edge: towards the vertex of triangle 0 not on the edge */
        int other=-1;
        for (int k=0;k<3;k++) if (tri[0][k]!=shared[0]&&tri[0][k]!=shared[1]) other=tri[0][k];
        double n=0;
        for (int k=0;k<3;k++) {dir[k]=sta[other].pos[k]-mid[k]; n+=dir[k]*dir[k];}
        for (int k=0;k<3;k++) dir[k]/=sqrt(n);
    }
    double pa[3],pb[3],da[3],db[3];
    nrtk_datum_info_t ia,ib;
    for (int k=0;k<3;k++) {pa[k]=mid[k]+1.0*dir[k]; pb[k]=mid[k]-1.0*dir[k];}   /* 1 m either side */
    EXPECT_EQ(nrtk_datum_shift(d,pa,t,da,&ia),NRTK_DSHIFT_TRI);
    EXPECT_EQ(nrtk_datum_shift(d,pb,t,db,&ib),NRTK_DSHIFT_TRI);
    /* different triangles ... */
    bool same=true;
    for (int k=0;k<3;k++) if (ia.pt[k]!=ib.pt[k]) same=false;
    EXPECT_FALSE(same);
    /* ... but the same value on the edge: 2 m apart in a field with gradients of
       1e-5 per metre differ by far less than 0.1 mm */
    for (int k=0;k<3;k++) EXPECT_NEAR(da[k],db[k],1E-4);
    nrtk_datum_free(d);
}

TEST(Datum, ExtrapolationNearestAndModes)
{
    nrtk_station_t sta[4]={sta_at("A",24.95,121.25,100),sta_at("B",24.70,120.92,50),
                           sta_at("C",24.62,121.30,600),sta_at("D",25.05,120.85,30)};
    LinField lf={24.85,121.05,{0.45,-0.30,-0.02},{0.004,0.002,0.0005},{-0.003,0.005,-0.0008}};
    for (auto &s:sta) lf.apply(s,0);
    nrtk_datum_def_t dd=def("TWD97_2010");
    nrtk_datum_t *d=nrtk_datum_create(&dd,0,sta,4,nullptr,0,NRTK_DATUM_INTERP_BARY,20.0,10.0);
    gtime_t t=t2026();
    /* 8 km east of the network: plane extended, still exact for a linear field */
    double llh[3]={24.80*D2R,121.40*D2R,100},pos[3],delta[3],enu[3],want[3];
    nrtk_datum_info_t info;
    pos2ecef(llh,pos);
    EXPECT_EQ(nrtk_datum_shift(d,pos,t,delta,&info),NRTK_DSHIFT_EXTRAP);
    EXPECT_GT(info.extrap_km,1.0);
    EXPECT_LT(info.extrap_km,20.0);
    lf.at(pos,want);
    for (int k=0;k<3;k++) EXPECT_NEAR(delta[k],want[k],2E-5);
    enu_of(pos,delta,enu);
    EXPECT_NEAR(info.enu[0],enu[0],1E-9);
    /* 60 km away: beyond the limit, nearest point's shift */
    llh[1]=121.95*D2R; pos2ecef(llh,pos);
    EXPECT_EQ(nrtk_datum_shift(d,pos,t,delta,&info),NRTK_DSHIFT_NEAREST);
    EXPECT_GT(info.extrap_km,20.0);
    EXPECT_GE(info.pt[0],0);
    EXPECT_EQ(info.pt[1],-1);
    double dn[3];
    nrtk_datum_point(d,info.pt[0],t,nullptr,dn);
    for (int k=0;k<3;k++) EXPECT_NEAR(delta[k],dn[k],1E-9);
    nrtk_datum_free(d);
    /* interp nearest: constant per Voronoi cell, even inside */
    d=nrtk_datum_create(&dd,0,sta,4,nullptr,0,NRTK_DATUM_INTERP_NEAREST,20.0,10.0);
    llh[0]=24.93*D2R; llh[1]=121.22*D2R; pos2ecef(llh,pos);
    EXPECT_EQ(nrtk_datum_shift(d,pos,t,delta,&info),NRTK_DSHIFT_NEAREST);
    EXPECT_EQ(info.pt[0],0);
    nrtk_datum_free(d);
    /* two points only: nearest; none: no shift */
    d=nrtk_datum_create(&dd,0,sta,2,nullptr,0,NRTK_DATUM_INTERP_BARY,20.0,10.0);
    EXPECT_EQ(nrtk_datum_ntri(d),0);
    EXPECT_EQ(nrtk_datum_shift(d,pos,t,delta,&info),NRTK_DSHIFT_NEAREST);
    nrtk_datum_free(d);
    d=nrtk_datum_create(&dd,1,sta,4,nullptr,0,NRTK_DATUM_INTERP_BARY,20.0,10.0);   /* datum 1: no coordinates */
    EXPECT_EQ(nrtk_datum_npoints(d),0);
    EXPECT_EQ(nrtk_datum_nmissing(d),4);
    EXPECT_EQ(nrtk_datum_shift(d,pos,t,delta,&info),NRTK_DSHIFT_NONE);
    for (int k=0;k<3;k++) EXPECT_EQ(delta[k],0.0);
    nrtk_datum_free(d);
}

TEST(Datum, InverseRoundTripAndTheoryFlag)
{
    nrtk_station_t sta[4]={sta_at("A",24.95,121.25,100),sta_at("B",24.70,120.92,50),
                           sta_at("C",24.62,121.30,600),sta_at("D",25.05,120.85,30)};
    double rnd[4][3]={{0.41,-0.33,0.02},{0.52,-0.21,-0.04},{0.38,-0.40,0.07},{0.60,-0.25,0.01}};
    for (int i=0;i<4;i++) set_official(sta[i],0,rnd[i][0],rnd[i][1],rnd[i][2],i==3?"theory":"control-tie");
    nrtk_datum_def_t dd=def("TWD97_2010");
    nrtk_datum_t *d=nrtk_datum_create(&dd,0,sta,4,nullptr,0,NRTK_DATUM_INTERP_BARY,20.0,10.0);
    gtime_t t=t2026();
    double llh[3]={24.86*D2R,121.12*D2R,150},pos[3],delta[3],pd[3],back[3];
    nrtk_datum_info_t info;
    pos2ecef(llh,pos);
    nrtk_datum_shift(d,pos,t,delta,&info);
    for (int k=0;k<3;k++) pd[k]=pos[k]+delta[k];
    EXPECT_EQ(nrtk_datum_inverse(d,pd,t,back),NRTK_DSHIFT_TRI);
    for (int k=0;k<3;k++) EXPECT_NEAR(back[k],pos[k],1E-3);
    /* theory-sourced point flagged only where it takes part */
    bool usesD=false;
    for (int k=0;k<3;k++) if (info.pt[k]==3&&info.w[k]>1E-9) usesD=true;
    EXPECT_EQ(info.theory?1:0,usesD?1:0);
    double dC[3];
    nrtk_datum_shift(d,sta[2].pos,t,dC,&info);  /* vertex C: weights of the others are 0 */
    EXPECT_FALSE(info.theory);
    nrtk_datum_free(d);
}

TEST(Datum, VelocityMovesTheShiftAndLeaveOneOut)
{
    nrtk_station_t sta[5]={sta_at("A",24.95,121.25,100),sta_at("B",24.70,120.92,50),
                           sta_at("C",24.62,121.30,600),sta_at("D",25.05,120.85,30),
                           sta_at("E",24.85,121.08,120)};
    LinField lf={24.85,121.05,{0.45,-0.30,-0.02},{0.004,0.002,0.0005},{-0.003,0.005,-0.0008}};
    for (auto &s:sta) lf.apply(s,0);
    nrtk_datum_def_t dd=def("TWD97_2010");
    nrtk_datum_t *d=nrtk_datum_create(&dd,0,sta,5,nullptr,0,NRTK_DATUM_INTERP_BARY,20.0,10.0);
    gtime_t t=t2026();
    /* leave-one-out of the interior point E in a linear field: zero residual */
    double li[3],de[3];
    nrtk_datum_info_t info;
    EXPECT_EQ(nrtk_datum_loo(d,4,t,li,&info),NRTK_DSHIFT_TRI);
    nrtk_datum_point(d,4,t,nullptr,de);
    for (int k=0;k<3;k++) EXPECT_NEAR(li[k],de[k],2E-5);
    for (int k=0;k<3;k++) EXPECT_NE(info.pt[k],4);      /* indices refer to the full set */
    nrtk_datum_free(d);
    /* station A moves 3 cm/yr in x since 2020.0: at 2026.78 its ITRF position is
       0.203 m further, the (static) official coordinates are not */
    sta[0].coord_epoch=2020.0; sta[0].vel[0]=0.03;
    d=nrtk_datum_create(&dd,0,sta,5,nullptr,0,NRTK_DATUM_INTERP_BARY,20.0,10.0);
    double dy=nrtk_decimal_year(t)-2020.0,pos[3],delta[3];
    nrtk_datum_point(d,0,t,pos,delta);
    EXPECT_NEAR(pos[0],sta[0].pos[0]+0.03*dy,1E-9);
    EXPECT_NEAR(delta[0],(sta[0].datum[0].pos[0]-sta[0].pos[0])-0.03*dy,1E-9);
    /* the field is no longer linear: E's leave-one-out picks up part of A's motion */
    EXPECT_EQ(nrtk_datum_loo(d,4,t,li,&info),NRTK_DSHIFT_TRI);
    nrtk_datum_point(d,4,t,nullptr,de);
    EXPECT_GT(fabs(li[0]-de[0]),0.01);
    nrtk_datum_free(d);
    /* fewer than 4 points: no leave-one-out */
    d=nrtk_datum_create(&dd,0,sta,3,nullptr,0,NRTK_DATUM_INTERP_BARY,20.0,10.0);
    EXPECT_EQ(nrtk_datum_loo(d,0,t,li,&info),NRTK_DSHIFT_NONE);
    nrtk_datum_free(d);
}

TEST(Datum, CalibrationPointsAndMissingStations)
{
    nrtk_station_t sta[3]={sta_at("A",24.95,121.25,100),sta_at("B",24.70,120.92,50),
                           sta_at("C",24.62,121.30,600)};
    set_official(sta[0],0,0.4,-0.3,0.0);
    set_official(sta[1],0,0.5,-0.2,0.0);
    /* C has no coordinates in the datum; a calibration point next to it has */
    nrtk_dpoint_t pt;
    memset(&pt,0,sizeof(pt));
    snprintf(pt.id,sizeof(pt.id),"GC01");
    nrtk_station_t tmp=sta_at("GC01",24.64,121.28,550);
    set_official(tmp,0,0.45,-0.35,0.0,"control-point");
    memcpy(pt.pos,tmp.pos,sizeof(pt.pos));
    pt.coord_epoch=2026.0;
    pt.datum[0]=tmp.datum[0];
    nrtk_datum_def_t dd=def("TWD97_2010");
    nrtk_datum_t *d=nrtk_datum_create(&dd,0,sta,3,&pt,1,NRTK_DATUM_INTERP_BARY,20.0,10.0);
    EXPECT_EQ(nrtk_datum_npoints(d),3);
    EXPECT_EQ(nrtk_datum_ntri(d),1);
    ASSERT_EQ(nrtk_datum_nmissing(d),1);
    EXPECT_EQ(nrtk_datum_missing(d,0),2);
    EXPECT_EQ(nrtk_datum_point_station(d,2),-1);
    EXPECT_STREQ(nrtk_datum_point_id(d,2),"GC01");
    EXPECT_STREQ(nrtk_datum_point_source(d,2),"control-point");
    /* a position inside A-B-GC01 is served (interpolated) although C is missing */
    double llh[3]={24.78*D2R,121.12*D2R,100},pos[3],delta[3];
    nrtk_datum_info_t info;
    pos2ecef(llh,pos);
    EXPECT_EQ(nrtk_datum_shift(d,pos,t2026(),delta,&info),NRTK_DSHIFT_TRI);
    /* GC01 is 3 km from C: C is covered, the network cell A-B-C is complete */
    int tri[3]={0,1,2};
    EXPECT_TRUE(nrtk_datum_missing_covered(d,0));
    EXPECT_EQ(nrtk_datum_nuncovered(d),0);
    EXPECT_FALSE(nrtk_datum_cell_incomplete(d,tri));
    nrtk_datum_free(d);
    /* with a 1 km coverage radius C is not covered: its cells are incomplete */
    d=nrtk_datum_create(&dd,0,sta,3,&pt,1,NRTK_DATUM_INTERP_BARY,20.0,1.0);
    EXPECT_FALSE(nrtk_datum_missing_covered(d,0));
    EXPECT_EQ(nrtk_datum_nuncovered(d),1);
    EXPECT_TRUE(nrtk_datum_cell_incomplete(d,tri));
    int tri2[3]={0,1,-1};
    EXPECT_FALSE(nrtk_datum_cell_incomplete(d,tri2));
    EXPECT_FALSE(nrtk_datum_cell_incomplete(d,nullptr));
    nrtk_datum_free(d);
}
