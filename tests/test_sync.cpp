/* unit tests: real-time epoch synchronization */
#include <gtest/gtest.h>
#include <thread>
#include "nrtk/realtime.hpp"

static std::vector<obsd_t> obs_at(gtime_t t, int sat)
{
    obsd_t o; memset(&o,0,sizeof(o)); o.time=t; o.sat=(uint8_t)sat;
    return {o};
}

TEST(EpochSync, CompleteEpochReleasedImmediately)
{
    nrtk::Clock clk;
    nrtk::EpochSync s(3,0.5,clk);
    gtime_t t=utc2gpst(timeget());
    for (int i=0;i<3;i++) s.push(i,obs_at(t,i+1));
    nrtk::SyncEpoch e;
    auto t0=std::chrono::steady_clock::now();
    ASSERT_TRUE(s.pop(e,1000));
    EXPECT_LT(std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count(),0.1);
    EXPECT_EQ(e.nsta,3);
    EXPECT_FALSE(e.raw[2].empty());
}

TEST(EpochSync, MissingStationAfterWindow)
{
    nrtk::Clock clk;
    nrtk::EpochSync s(3,0.3,clk);
    gtime_t t=utc2gpst(timeget());
    for (int i=0;i<3;i++) s.push(i,obs_at(timeadd(t,-1),i+1));   /* all alive */
    nrtk::SyncEpoch e;
    ASSERT_TRUE(s.pop(e,1000));
    s.push(0,obs_at(t,1)); s.push(1,obs_at(t,2));                /* station 2 late */
    auto t0=std::chrono::steady_clock::now();
    ASSERT_TRUE(s.pop(e,2000));
    double w=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    EXPECT_GT(w,0.25);
    EXPECT_LT(w,0.6);
    EXPECT_EQ(e.nsta,2);
    EXPECT_TRUE(e.raw[2].empty());
    s.push(2,obs_at(t,3));                                       /* too late: dropped */
    EXPECT_EQ(s.dropped_late(),1);
    EXPECT_FALSE(s.pop(e,100));
}

TEST(EpochSync, OrderedRelease)
{
    nrtk::Clock clk;
    nrtk::EpochSync s(2,0.2,clk);
    gtime_t t=utc2gpst(timeget());
    s.push(0,obs_at(timeadd(t,1),1));
    s.push(0,obs_at(t,1));
    s.push(1,obs_at(t,2));
    s.push(1,obs_at(timeadd(t,1),2));
    nrtk::SyncEpoch a,b;
    ASSERT_TRUE(s.pop(a,500));
    ASSERT_TRUE(s.pop(b,500));
    EXPECT_LT(timediff(a.t,b.t),0.0);
}
