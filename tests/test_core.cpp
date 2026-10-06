/* unit tests: troposphere, ionosphere split, interpolation, cells, antenna */
#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include "nrtk/common.h"
#include "nrtk/tropo.h"
#include "nrtk/netpos.h"
#include "nrtk/atmos.h"
#include "nrtk/cellmgr.h"
#include "nrtk/antenna.h"

#define SQR(x) ((x)*(x))

/* IERS Conventions test case of GMF (gmf.f header) */
TEST(Tropo, GmfIersTestCase)
{
    double gh,gw;
    nrtk_gmf(55055.0,0.6708665767,-1.393397187,844.715,1.278564131,&gh,&gw);
    EXPECT_NEAR(gh,3.425245519339138678,1E-12);
    EXPECT_NEAR(gw,3.449589116182419257,1E-12);
}

TEST(Tropo, SlantDelayPlausible)
{
    double pos[3]={24.8*D2R,121.0*D2R,100.0},azel[2]={0.0,90.0*D2R},zhd,zwd;
    gtime_t t=epoch2time(std::vector<double>{2026,9,20,4,0,0}.data());
    nrtk_trop_zenith(pos,0.7,&zhd,&zwd);
    EXPECT_NEAR(zhd,2.28,0.05);
    EXPECT_GT(zwd,0.05);
    double z=nrtk_trop_slant(t,pos,azel);
    EXPECT_NEAR(z,zhd+zwd,0.01);
    azel[1]=10.0*D2R;
    double s=nrtk_trop_slant(t,pos,azel);
    EXPECT_GT(s/z,5.0);           /* mapping at 10 deg ~5.6 */
    EXPECT_LT(s/z,6.2);
    pos[2]=600.0;                 /* height dependence */
    azel[1]=90.0*D2R;
    EXPECT_LT(nrtk_trop_slant(t,pos,azel),z-0.1);
}

/* residual split: da = D - (f1/fa)^2 I, db = D - (f1/fb)^2 I */
TEST(Netpos, SplitIonoRoundTrip)
{
    const double f[][2]={{FREQL1,FREQL2},{FREQL1,FREQE5b},{FREQ1_CMP,FREQ3_CMP},{FREQL2,FREQL5}};
    for (auto &p:f) {
        double I=0.437,D=-0.123,ga=SQR(FREQL1/p[0]),gb=SQR(FREQL1/p[1]),Ie,De;
        nrtk_split_iono(p[0],p[1],D-ga*I,D-gb*I,&Ie,&De);
        EXPECT_NEAR(Ie,I,1E-12);
        EXPECT_NEAR(De,D,1E-12);
    }
    /* applying the split back to L1/L2 reproduces the residuals exactly */
    double d1=0.031,d2=-0.012,I,D;
    nrtk_split_iono(FREQL1,FREQL2,d1,d2,&I,&D);
    EXPECT_NEAR(D-I,d1,1E-12);
    EXPECT_NEAR(D-SQR(FREQL1/FREQL2)*I,d2,1E-12);
}

TEST(Atmos, PlaneFitLimExact)
{
    double n[2]={10.0,-20.0},e[2]={30.0,5.0},y[2],a,b,rms;
    for (int i=0;i<2;i++) y[i]=0.002*n[i]-0.001*e[i];
    ASSERT_TRUE(nrtk_plane_fit(n,e,y,2,&a,&b,&rms));
    EXPECT_NEAR(a,0.002,1E-12);
    EXPECT_NEAR(b,-0.001,1E-12);
    EXPECT_EQ(rms,0.0);
}

TEST(Atmos, PlaneFitLsmResidual)
{
    double n[4]={10,-20,5,-8},e[4]={30,5,-25,-12},y[4],a,b,rms;
    for (int i=0;i<4;i++) y[i]=0.001*n[i]+0.0005*e[i];
    y[2]+=0.01;                       /* non-planar part */
    ASSERT_TRUE(nrtk_plane_fit(n,e,y,4,&a,&b,&rms));
    EXPECT_GT(rms,0.003);
    EXPECT_NEAR(a,0.001,5E-4);
}

TEST(Atmos, PlaneFitRejectsCollinear)
{
    double n[3]={1,2,3},e[3]={2,4,6},y[3]={0.1,0.2,0.3},a,b,rms;
    EXPECT_FALSE(nrtk_plane_fit(n,e,y,3,&a,&b,&rms));
    EXPECT_FALSE(nrtk_plane_fit(n,e,y,1,&a,&b,&rms));
}

/* Delaunay: empty circumcircle property on random points */
TEST(Cellmgr, DelaunayEmptyCircumcircle)
{
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-50,50);
    double x[12],y[12];
    for (int i=0;i<12;i++) {x[i]=u(rng); y[i]=u(rng);}
    int tri[64][3];
    int nt=nrtk_delaunay(x,y,12,tri,64);
    ASSERT_GE(nt,10);
    for (int k=0;k<nt;k++) {
        double ax=x[tri[k][0]],ay=y[tri[k][0]],bx=x[tri[k][1]],by=y[tri[k][1]],cx=x[tri[k][2]],cy=y[tri[k][2]];
        double d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by));
        double ux=((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/d;
        double uy=((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/d;
        double r2=(ax-ux)*(ax-ux)+(ay-uy)*(ay-uy);
        for (int i=0;i<12;i++) {
            if (i==tri[k][0]||i==tri[k][1]||i==tri[k][2]) continue;
            EXPECT_GT((x[i]-ux)*(x[i]-ux)+(y[i]-uy)*(y[i]-uy),r2*(1-1E-9));
        }
    }
    /* square: two triangles */
    double sx[4]={0,10,10,0},sy[4]={0,0,10,10};
    EXPECT_EQ(nrtk_delaunay(sx,sy,4,tri,64),2);
}

static nrtk_station_t make_sta(const char *id, double lat, double lon, double h)
{
    nrtk_station_t s;
    memset(&s,0,sizeof(s));
    snprintf(s.id,sizeof(s.id),"%s",id);
    double p[3]={lat*D2R,lon*D2R,h};
    pos2ecef(p,s.pos);
    s.enabled=1;
    return s;
}

TEST(Cellmgr, SelectMasterAuxExtrapReject)
{
    nrtk_station_t sta[4]={make_sta("A",24.95,121.25,100),make_sta("B",24.70,120.92,50),
                           make_sta("C",24.62,121.30,600),make_sta("D",25.05,120.85,30)};
    nrtk_cellmgr_opt_t opt;
    nrtk_cellmgr_opt_default(&opt);
    nrtk_cellmgr_t *cm=nrtk_cellmgr_create(&opt,sta,4);
    nrtk_cell_t cell;
    double p[3],v[3];
    /* inside A-B-C, closest to A */
    p[0]=24.85*D2R; p[1]=121.15*D2R; p[2]=100; pos2ecef(p,v);
    ASSERT_TRUE(nrtk_cellmgr_select(cm,v,&cell));
    EXPECT_EQ(cell.master,0);
    EXPECT_FALSE(cell.extrap);
    EXPECT_GE(cell.naux,2);
    /* outside the network (east), still within 70 km of A */
    p[0]=24.90*D2R; p[1]=121.55*D2R; pos2ecef(p,v);
    ASSERT_TRUE(nrtk_cellmgr_select(cm,v,&cell));
    EXPECT_TRUE(cell.extrap);
    /* far away: rejected */
    p[0]=23.5*D2R; p[1]=121.5*D2R; pos2ecef(p,v);
    EXPECT_FALSE(nrtk_cellmgr_select(cm,v,&cell));
    EXPECT_TRUE(cell.reject);
    /* unhealthy master: another station becomes master */
    p[0]=24.85*D2R; p[1]=121.15*D2R; pos2ecef(p,v);
    EXPECT_TRUE(nrtk_cellmgr_set_health(cm,0,0));
    ASSERT_TRUE(nrtk_cellmgr_select(cm,v,&cell));
    EXPECT_NE(cell.master,0);
    /* fewer than 3 healthy stations: no cell */
    nrtk_cellmgr_set_health(cm,1,0);
    EXPECT_FALSE(nrtk_cellmgr_select(cm,v,&cell));
    nrtk_cellmgr_free(cm);
}

TEST(Common, SignalHelpers)
{
    EXPECT_EQ(nrtk_sysmask("GPS,GAL"),(1<<NRTK_GPS)|(1<<NRTK_GAL));
    EXPECT_EQ(nrtk_sysmask("[BDS, QZSS]"),(1<<NRTK_BDS)|(1<<NRTK_QZS));
    EXPECT_TRUE(nrtk_sig_in_list(obs2code("2W"),"1C 2W 5Q"));
    EXPECT_FALSE(nrtk_sig_in_list(obs2code("2L"),"1C 2W 5Q"));
    EXPECT_DOUBLE_EQ(nrtk_code2freq(SYS_GAL,obs2code("7Q")),FREQE5b);
    EXPECT_DOUBLE_EQ(nrtk_code2freq(SYS_CMP,obs2code("6I")),FREQ3_CMP);
    nrtk_sigplan_t plan[NRTK_NSYS];
    nrtk_sigplan_default(plan);
    EXPECT_DOUBLE_EQ(plan[NRTK_GAL].pfreq[1],FREQE5b);   /* E1/E5b processing pair */
    EXPECT_TRUE(nrtk_sig_in_list(obs2code("6I"),plan[NRTK_BDS].proc[1]));
}

/* minimal ANTEX with one receiver antenna (two frequencies) and one satellite */
static std::string write_atx()
{
    std::string fn=testing::TempDir()+"nrtk_test.atx";
    std::ofstream f(fn);
    auto line=[&](const std::string &c, const std::string &lab){
        char b[128]; snprintf(b,sizeof(b),"%-60s%-20s\n",c.c_str(),lab.c_str()); f<<b;};
    line("     1.4            M","ANTEX VERSION / SYST");
    line("","END OF HEADER");
    line("","START OF ANTENNA");
    line("BLOCK IIF           G01                 G063      2010-022A","TYPE / SERIAL NO");
    line("     0.0  17.0   1.0","ZEN1 / ZEN2 / DZEN");
    line("   G01","START OF FREQUENCY");
    line("      0.00      0.00   1000.00","NORTH / EAST / UP");
    line("   G01","END OF FREQUENCY");
    line("","END OF ANTENNA");
    line("","START OF ANTENNA");
    line("TESTANT         NONE","TYPE / SERIAL NO");
    line("     0.0","DAZI");
    line("     0.0  90.0  10.0","ZEN1 / ZEN2 / DZEN");
    line("   G01","START OF FREQUENCY");
    line("      1.00      2.00     60.00","NORTH / EAST / UP");
    f<<"   NOAZI    0.00    1.00    2.00    3.00    4.00    5.00    6.00    7.00    8.00    9.00\n";
    line("   G01","END OF FREQUENCY");
    line("   E05","START OF FREQUENCY");
    line("      0.00      0.00     50.00","NORTH / EAST / UP");
    f<<"   NOAZI    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\n";
    line("   E05","END OF FREQUENCY");
    line("","END OF ANTENNA");
    return fn;
}

TEST(Antenna, ParseAndCorrect)
{
    nrtk_atx_t atx={0};
    std::string fn=write_atx();
    ASSERT_TRUE(nrtk_atx_read(fn.c_str(),&atx));
    EXPECT_EQ(atx.n,1);                               /* satellite antenna skipped */
    const nrtk_ant_t *a=nrtk_atx_find(&atx,"TESTANT NONE");
    ASSERT_NE(a,nullptr);
    EXPECT_EQ(nrtk_atx_find(&atx,"TESTANT"),a);       /* radome defaults to NONE */
    EXPECT_EQ(nrtk_atx_find(&atx,"TESTANT LEIT"),a);  /* falls back to NONE */
    EXPECT_EQ(nrtk_atx_find(&atx,"OTHER NONE"),nullptr);
    double azel[2]={0.0,90.0*D2R},del[3]={0,0,0.1};
    /* zenith: -(up+ARP) + PCV(z=0) */
    EXPECT_NEAR(nrtk_ant_dr(a,SYS_GPS,FREQL1,azel,del,1),-(0.060+0.1),1E-9);
    /* elevation 45 deg -> zenith 45: PCV interpolated 4.5 mm; north/east projections */
    azel[0]=0.0; azel[1]=45*D2R;
    double e=-(0.001*cos(45*D2R)+(0.060+0.1)*sin(45*D2R))+0.0045;
    EXPECT_NEAR(nrtk_ant_dr(a,SYS_GPS,FREQL1,azel,del,1),e,1E-9);
    EXPECT_NEAR(nrtk_ant_dr(a,SYS_GPS,FREQL1,azel,del,0),e-0.0045,1E-9);
    /* Galileo E5a uses E05; BDS B1I (no Cxx) falls back to G01 */
    azel[1]=90*D2R;
    EXPECT_NEAR(nrtk_ant_dr(a,SYS_GAL,FREQL5,azel,del,1),-(0.050+0.1),1E-9);
    EXPECT_NEAR(nrtk_ant_dr(a,SYS_CMP,FREQ1_CMP,azel,del,1),-(0.060+0.1),1E-9);
    /* no antenna: eccentricity only */
    EXPECT_NEAR(nrtk_ant_dr(nullptr,SYS_GPS,FREQL1,azel,del,1),-0.1,1E-12);
    nrtk_atx_free(&atx);
}
