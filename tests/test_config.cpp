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
