/* unit tests: YAML configuration */
#include <gtest/gtest.h>
#include <fstream>
#include "nrtk/config.h"

static std::string write(const std::string &name, const std::string &text)
{
    std::string fn=testing::TempDir()+name;
    std::ofstream(fn)<<text;
    return fn;
}
static const char *base_yaml=R"(
network:
  frame: ITRF2020
  epoch: 2026.75
  max_baseline_km: 70
stations:
  - id: STA1
    ecef: [-2956000.123, 5000000.456, 2600000.789]
    antenna: SEPCHOKE_B3E6  NONE
    arp_height_m: 0.0
    stream: ntrip://user:pass@caster:2101/STA1
  - id: STA2
    llh: [24.8, 121.1, 120.0]
    antenna: "LEIAR25.R4 LEIT"
    arp_height_m: 0.105
solver:
  systems: [GPS, GAL, BDS, QZS]
  elevation_mask_deg: 10
  ar_ratio: 3.0
  iono_prior_ppm: 0.5
  engine: rtklib
output:
  msm: 7
  rate_hz: 1
  vrs_rebuild_km: 5
)";

TEST(Config, SpecExampleLoads)
{
    nrtk_config_t *c=new nrtk_config_t;
    char err[512]="";
    nrtk_config_default(c);
    ASSERT_TRUE(nrtk_config_load(write("c1.yaml",base_yaml).c_str(),c,err,sizeof(err)))<<err;
    EXPECT_EQ(c->nsta,2);
    EXPECT_STREQ(c->sta[0].id,"STA1");
    EXPECT_STREQ(c->sta[0].antenna,"SEPCHOKE_B3E6  NONE");
    EXPECT_STREQ(c->sta[0].source,"ntrip://user:pass@caster:2101/STA1");
    EXPECT_NEAR(c->sta[1].antdel[2],0.105,1E-12);
    double pos[3]; ecef2pos(c->sta[1].pos,pos);
    EXPECT_NEAR(pos[0]*R2D,24.8,1E-9);
    EXPECT_NEAR(pos[2],120.0,1E-4);
    EXPECT_EQ(c->netpos.sysmask,0xF);
    EXPECT_NEAR(c->netpos.elmin*R2D,10.0,1E-9);
    EXPECT_EQ(c->netpos.engine,NRTK_ENGINE_RTKLIB);
    EXPECT_EQ(c->netpos.closure,1);             /* default on */
    EXPECT_EQ(c->msm,7);
    EXPECT_NEAR(c->vrs_rebuild_km,5.0,1E-12);
    EXPECT_NEAR(c->cell.max_baseline_km,70.0,1E-12);
    EXPECT_EQ(nrtk_config_validate(c,err,sizeof(err)),0)<<err;
    delete c;
}

TEST(Config, ClosureCheckSwitch)
{
    nrtk_config_t *c=new nrtk_config_t;
    char err[512]="";
    nrtk_config_default(c);
    std::string y=std::string(base_yaml)+"";
    y.replace(y.find("  engine: rtklib\n"),16,"  engine: native\n  closure_check: false\n");
    ASSERT_TRUE(nrtk_config_load(write("c_closure.yaml",y).c_str(),c,err,sizeof(err)))<<err;
    EXPECT_EQ(c->netpos.engine,NRTK_ENGINE_NATIVE);
    EXPECT_EQ(c->netpos.closure,0);
    delete c;
}

TEST(Config, ValidationErrors)
{
    std::string y=base_yaml;
    y.replace(y.find("msm: 7"),6,"msm: 5");
    y.replace(y.find("id: STA2"),8,"id: STA1");
    y+="offline:\n  vrs:\n    - {name: x, leave_out: NOPE}\n";
    nrtk_config_t *c=new nrtk_config_t;
    char err[1024]="";
    nrtk_config_default(c);
    ASSERT_TRUE(nrtk_config_load(write("c2.yaml",y).c_str(),c,err,sizeof(err)))<<err;
    EXPECT_EQ(nrtk_config_validate(c,err,sizeof(err)),3);
    EXPECT_NE(std::string(err).find("msm"),std::string::npos);
    EXPECT_NE(std::string(err).find("duplicate"),std::string::npos);
    EXPECT_NE(std::string(err).find("NOPE"),std::string::npos);
    delete c;
}

TEST(Config, BadInputs)
{
    nrtk_config_t *c=new nrtk_config_t;
    char err[512]="";
    nrtk_config_default(c);
    EXPECT_FALSE(nrtk_config_load("/nonexistent/x.yaml",c,err,sizeof(err)));
    EXPECT_FALSE(nrtk_config_load(write("c3.yaml","stations:\n  - {ecef: [1,2,3]}\n").c_str(),c,err,sizeof(err)));
    EXPECT_FALSE(nrtk_config_load(write("c4.yaml","solver: {engine: magic}\n").c_str(),c,err,sizeof(err)));
    gtime_t t;
    EXPECT_TRUE(nrtk_parse_time("2026-09-20 04:00:00",&t));
    EXPECT_FALSE(nrtk_parse_time("yesterday",&t));
    delete c;
}

/* output datums: definitions, station coordinates, calibration points, mounts */
#include "nrtk/datum.h"
TEST(Config, OutputDatums)
{
    std::string pts=write("pts.yaml",
        "points:\n"
        "  - {id: GC01, llh: [24.75, 121.00, 80.0], epoch: 2026.5, velocity: [0.01, 0.02, 0.0],\n"
        "     datum: {TWD97_2010: {llh: [24.750003, 121.000004, 79.98], source: control-point}}}\n");
    std::string y=std::string(base_yaml);
    y.replace(y.find("  max_baseline_km: 70\n"),22,
        "  max_baseline_km: 70\n"
        "  output_datums:\n"
        "    - {name: TWD97_2010, label: \"TWD97[2010]\", version: cal-2026.10}\n"
        "    - {name: TWD97_2020}\n"
        "  datum_interp: barycentric\n"
        "  datum_extrap_max_km: 15\n"
        "  datum_points_file: "+pts+"\n");
    y.replace(y.find("    arp_height_m: 0.0\n"),22,
        "    arp_height_m: 0.0\n"
        "    datum:\n"
        "      TWD97_2010: {ecef: [-2956000.600, 5000000.200, 2600000.500], source: control-tie, observed: 2026-09-20, sigma_m: 0.015}\n"
        "      TWD97_2020: {ecef: [-2956000.400, 5000000.300, 2600000.600]}\n");
    y+="server:\n  datum_mounts:\n    - {datum: TWD97_2010}\n    - {datum: TWD97_2020, msm7: VRS_2020, nearest: \"\"}\n";
    nrtk_config_t *c=new nrtk_config_t;
    char err[1024]="";
    nrtk_config_default(c);
    ASSERT_TRUE(nrtk_config_load(write("c_datum.yaml",y).c_str(),c,err,sizeof(err)))<<err;
    ASSERT_EQ(c->ndatum,2);
    EXPECT_STREQ(c->datum[0].name,"TWD97_2010");
    EXPECT_STREQ(c->datum[0].label,"TWD97[2010]");
    EXPECT_STREQ(c->datum[0].version,"cal-2026.10");
    EXPECT_STREQ(c->datum[0].frame,"ITRF94");        /* theory defaults by name */
    EXPECT_DOUBLE_EQ(c->datum[0].epoch,2010.0);
    EXPECT_STREQ(c->datum[1].label,"TWD97_2020");
    EXPECT_STREQ(c->datum[1].frame,"ITRF2014");
    EXPECT_DOUBLE_EQ(c->datum[1].epoch,2020.0);
    EXPECT_EQ(nrtk_config_datum(c,"TWD97_2020"),1);
    EXPECT_EQ(nrtk_config_datum(c,"nope"),-1);
    EXPECT_EQ(c->datum_interp,NRTK_DATUM_INTERP_BARY);
    EXPECT_DOUBLE_EQ(c->datum_extrap_max_km,15.0);
    /* station STA1: both datums; STA2: none */
    EXPECT_TRUE(c->sta[0].datum[0].valid);
    EXPECT_TRUE(c->sta[0].datum[1].valid);
    EXPECT_FALSE(c->sta[1].datum[0].valid);
    EXPECT_NEAR(c->sta[0].datum[0].pos[0],-2956000.600,1E-9);
    EXPECT_NEAR(c->sta[0].datum[0].sigma,0.015,1E-12);
    EXPECT_STREQ(c->sta[0].datum[0].source,"control-tie");
    EXPECT_STREQ(c->sta[0].datum[0].observed,"2026-09-20");
    /* calibration point from the file */
    ASSERT_EQ(c->ndpoint,1);
    EXPECT_STREQ(c->dpoint[0].id,"GC01");
    EXPECT_DOUBLE_EQ(c->dpoint[0].coord_epoch,2026.5);
    EXPECT_DOUBLE_EQ(c->dpoint[0].vel[1],0.02);
    EXPECT_TRUE(c->dpoint[0].datum[0].valid);
    EXPECT_FALSE(c->dpoint[0].datum[1].valid);
    EXPECT_LT(nrtk_dist(c->dpoint[0].datum[0].pos,c->dpoint[0].pos),2.0);
    /* mounts: default names <mount>_<datum>, explicit overrides, empty disables */
    ASSERT_EQ(c->ndmount,2);
    EXPECT_EQ(c->dmount[0].datum,0);
    EXPECT_STREQ(c->dmount[0].msm7,"VRS_MSM7_TWD97_2010");
    EXPECT_STREQ(c->dmount[0].msm4,"VRS_MSM4_TWD97_2010");
    EXPECT_STREQ(c->dmount[0].near,"NEAR_TWD97_2010");
    EXPECT_STREQ(c->dmount[1].msm7,"VRS_2020");
    EXPECT_STREQ(c->dmount[1].msm4,"VRS_MSM4_TWD97_2020");
    EXPECT_STREQ(c->dmount[1].near,"");
    EXPECT_EQ(nrtk_config_validate(c,err,sizeof(err)),0)<<err;
    /* shift field of datum 0: STA1 + GC01 (2 points), STA2 missing */
    nrtk_datum_t *d=nrtk_datum_from_config(c,0);
    ASSERT_NE(d,nullptr);
    EXPECT_EQ(nrtk_datum_npoints(d),2);
    EXPECT_EQ(nrtk_datum_nmissing(d),1);
    EXPECT_EQ(nrtk_datum_missing(d,0),1);
    nrtk_datum_free(d);
    /* validation: official coordinates metres away from ITRF (wrong version/unit) */
    c->sta[0].datum[1].pos[2]+=7.0;
    EXPECT_EQ(nrtk_config_validate(c,err,sizeof(err)),1);
    EXPECT_NE(std::string(err).find("TWD97_2020"),std::string::npos);
    delete c;
    /* load errors: unknown datum under a station, unknown datum in a mount, bad name */
    c=new nrtk_config_t;
    nrtk_config_default(c);
    std::string y2=std::string(base_yaml);
    y2.replace(y2.find("    arp_height_m: 0.0\n"),22,"    arp_height_m: 0.0\n    datum: {TWD97_2010: {ecef: [1,2,3]}}\n");
    EXPECT_FALSE(nrtk_config_load(write("c_d2.yaml",y2).c_str(),c,err,sizeof(err)));
    EXPECT_NE(std::string(err).find("output_datums"),std::string::npos);
    std::string y3=std::string(base_yaml)+"server:\n  datum_mounts: [{datum: X}]\n";
    EXPECT_FALSE(nrtk_config_load(write("c_d3.yaml",y3).c_str(),c,err,sizeof(err)));
    std::string y4=std::string(base_yaml);
    y4.replace(y4.find("  max_baseline_km: 70\n"),22,"  max_baseline_km: 70\n  output_datums: [{name: \"TWD97[2010]\"}]\n");
    EXPECT_FALSE(nrtk_config_load(write("c_d4.yaml",y4).c_str(),c,err,sizeof(err)));
    delete c;
}
