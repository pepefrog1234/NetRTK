/* unit tests: RTCM 3 output (1006, 1033, MSM7/MSM4 round trip with splitting) */
#include <gtest/gtest.h>
#include <cmath>
#include <map>
#include <vector>
#include "nrtk/rtcm_out.h"

static std::vector<rtcm_t *> decode_all(const uint8_t *buf, int n, gtime_t t,
                                        std::vector<std::vector<obsd_t>> *epochs,
                                        int *nmsg=nullptr)
{
    rtcm_t *r=new rtcm_t;
    init_rtcm(r);
    r->time=t;
    int m=0;
    for (int i=0;i<n;i++) {
        int ret=input_rtcm3(r,buf[i]);
        if (r->nbyte==0&&r->len>0) {}
        if (ret==1&&epochs) epochs->push_back(std::vector<obsd_t>(r->obs.data,r->obs.data+r->obs.n));
        if (ret!=0) m++;
    }
    if (nmsg) *nmsg=m;
    return {r};
}

TEST(Rtcm, Encode1006RoundTrip)
{
    double pos[3]={-3001875.5703,4946937.8364,2674117.0167};
    uint8_t buf[64];
    int n=nrtk_rtcm_encode_1006(3001,pos,0.0,1,0,0xF,buf,sizeof(buf));
    ASSERT_EQ(n,3+21+3);                     /* 168-bit body (DF002..DF028) */
    EXPECT_EQ(getbitu(buf,24,12),1006u);
    EXPECT_EQ(getbitu(buf,36,12),3001u);
    EXPECT_EQ(getbitu(buf,24+12+12+6+3,1),1u);   /* DF141 non-physical */
    rtcm_t *r=new rtcm_t; init_rtcm(r);
    int ret=0;
    for (int i=0;i<n;i++) ret=input_rtcm3(r,buf[i]);
    EXPECT_EQ(ret,5);
    EXPECT_EQ(r->staid,3001);
    for (int i=0;i<3;i++) EXPECT_NEAR(r->sta.pos[i],pos[i],1E-4);
    free_rtcm(r); delete r;
}

TEST(Rtcm, Station1033AntennaDescriptor)
{
    nrtk_rtcmout_t o;
    double pos[3]={-2997379.911,4959023.866,2657073.102};
    ASSERT_TRUE(nrtk_rtcmout_init(&o,1234,7,pos,1));
    uint8_t buf[512];
    int n=nrtk_rtcmout_station(&o,buf,sizeof(buf));
    rtcm_t *r=new rtcm_t; init_rtcm(r);
    for (int i=0;i<n;i++) input_rtcm3(r,buf[i]);
    EXPECT_STREQ(r->sta.antdes,"ADVNULLANTENNA");
    EXPECT_EQ(r->staid,1234);
    free_rtcm(r); delete r;
    nrtk_rtcmout_free(&o);
}

/* synthetic epoch: GPS 10 sats x 3 signals, BDS 18 sats x 4 signals (72 cells) */
static std::vector<obsd_t> make_epoch(gtime_t t)
{
    std::vector<obsd_t> v;
    struct {int sys; int nsat; int prn0; std::vector<const char *> sig;} grp[]={
        {SYS_GPS,10,1,{"1C","2W","5Q"}},{SYS_CMP,18,19,{"2I","6I","1P","5P"}},
        {SYS_GAL,8,1,{"1C","7Q","5Q"}}};
    for (auto &g:grp) for (int k=0;k<g.nsat;k++) {
        obsd_t o; memset(&o,0,sizeof(o));
        o.time=t; o.sat=(uint8_t)satno(g.sys,g.prn0+k); o.rcv=1;
        double rho=2.05E7+k*123456.789+(g.sys==SYS_CMP?1.5E6:0.0);
        for (size_t j=0;j<g.sig.size();j++) {
            o.code[j]=obs2code(g.sig[j]);
            double f=code2freq(g.sys,o.code[j],0),lam=CLIGHT/f;
            o.P[j]=rho+0.37*j;
            o.L[j]=(rho-0.21*j)/lam+12345.0+0.3*k;
            o.D[j]=(float)(-(234.5+k*10.0)/lam);    /* one range rate per satellite */
            o.SNR[j]=(float)(40.0+j);
        }
        v.push_back(o);
    }
    return v;
}

TEST(Rtcm, MsmRoundTripWithSplitting)
{
    for (int msm: {7,4}) {
        double pos[3]={-2997379.911,4959023.866,2657073.102};
        double ep[6]={2026,9,20,4,0,0};
        gtime_t t=epoch2time(ep);
        nrtk_rtcmout_t o;
        ASSERT_TRUE(nrtk_rtcmout_init(&o,77,msm,pos,1));
        std::vector<obsd_t> in=make_epoch(t);
        std::vector<uint8_t> buf(65536);
        int n=nrtk_rtcmout_obs(&o,t,in.data(),(int)in.size(),buf.data(),(int)buf.size());
        ASSERT_GT(n,0);
        /* count messages and check the multiple-message bit (only last is 0) */
        std::vector<int> types,sync;
        for (int i=0;i<n;) {
            ASSERT_EQ(buf[i],0xD3);
            int len=getbitu(buf.data()+i,14,10);
            types.push_back(getbitu(buf.data()+i,24,12));
            sync.push_back(getbitu(buf.data()+i,24+12+12+30,1));
            i+=len+6;
        }
        int nbds=0;
        for (int ty:types) if (ty==1120+msm) nbds++;
        EXPECT_EQ(nbds,2);                        /* 72 BDS cells split into 2 messages */
        for (size_t i=0;i+1<sync.size();i++) EXPECT_EQ(sync[i],1);
        EXPECT_EQ(sync.back(),0);
        /* decode with RTKLIB and compare */
        std::vector<std::vector<obsd_t>> eps;
        decode_all(buf.data(),n,t,&eps);
        ASSERT_EQ(eps.size(),1u);
        std::map<int,obsd_t> dec;
        for (auto &d:eps[0]) dec[d.sat]=d;
        ASSERT_EQ(dec.size(),in.size());
        for (auto &s:in) {
            const obsd_t &d=dec[s.sat];
            for (int j=0;j<NFREQ+NEXOBS;j++) {
                if (!s.code[j]) continue;
                int k=-1;
                for (int m=0;m<NFREQ+NEXOBS;m++) if (d.code[m]==s.code[j]) k=m;
                ASSERT_GE(k,0);
                double tolP=msm==7?0.002:0.02,tolL=msm==7?0.002:0.02;
                EXPECT_NEAR(d.P[k],s.P[j],tolP);
                double dl=d.L[k]-s.L[j];
                EXPECT_NEAR(dl-floor(dl+0.5),0.0,tolL);   /* integer offset only */
                char id[8]; satno2id(s.sat,id);
                if (msm==7) EXPECT_NEAR(d.D[k],s.D[j],0.01)<<id<<" "<<code2obs(s.code[j])<<" msm"<<msm;
                EXPECT_NEAR(d.SNR[k],s.SNR[j],msm==7?0.07:1.0)<<id<<" "<<code2obs(s.code[j])<<" msm"<<msm;
            }
        }
        nrtk_rtcmout_free(&o);
    }
}

TEST(Rtcm, LockTimeResetOnSlipFlag)
{
    double pos[3]={-2997379.911,4959023.866,2657073.102},ep[6]={2026,9,20,4,0,0};
    gtime_t t=epoch2time(ep);
    nrtk_rtcmout_t o;
    nrtk_rtcmout_init(&o,5,7,pos,1);
    rtcm_t *r=new rtcm_t; init_rtcm(r); r->time=t;
    std::vector<uint8_t> buf(65536);
    std::vector<obsd_t> in;
    int slipped=0,cont=0;
    for (int k=0;k<20;k++) {
        gtime_t tk=timeadd(t,k);
        in=make_epoch(tk);
        if (k==15) in[0].LLI[0]=LLI_SLIP;       /* reset lock of G01 L1 */
        int n=nrtk_rtcmout_obs(&o,tk,in.data(),(int)in.size(),buf.data(),(int)buf.size());
        for (int i=0;i<n;i++) {
            if (input_rtcm3(r,buf[i])==1&&k==15) {
                for (int m=0;m<r->obs.n;m++) if (r->obs.data[m].sat==in[0].sat) {
                    slipped=r->obs.data[m].LLI[0]&LLI_SLIP;
                    cont=r->obs.data[m].LLI[1]&LLI_SLIP;
                }
            }
        }
    }
    EXPECT_TRUE(slipped);        /* decoder sees the lock time drop */
    EXPECT_FALSE(cont);          /* other signals continuous */
    free_rtcm(r); delete r;
    nrtk_rtcmout_free(&o);
}

/* DF025..DF027 are 0.1 mm integers: the encoder rounds (no truncation bias) */
TEST(Rtcm, Encode1006RoundsPosition)
{
    double pos[3]={-3001875.57036,4946937.83644,2674117.01667};   /* x.xxxx6, x.xxxx4, x.xxxx7 */
    uint8_t buf[64];
    int n=nrtk_rtcm_encode_1006(1,pos,0.0,1,0,0xF,buf,sizeof(buf));
    ASSERT_GT(n,0);
    rtcm_t *r=new rtcm_t; init_rtcm(r);
    int ret=0;
    for (int i=0;i<n;i++) ret=input_rtcm3(r,buf[i]);
    EXPECT_EQ(ret,5);
    EXPECT_NEAR(r->sta.pos[0],-3001875.5704,1E-6);
    EXPECT_NEAR(r->sta.pos[1], 4946937.8364,1E-6);
    EXPECT_NEAR(r->sta.pos[2], 2674117.0167,1E-6);
    free_rtcm(r); delete r;
}
