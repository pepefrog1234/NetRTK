/*------------------------------------------------------------------------------
 * test_ntrip_front.cpp : NTRIP front end tests (GGA parser, base64, caster)
 *
 * Server tests run on 127.0.0.1 port 0 with raw client sockets, a fake
 * StreamSession that echoes the published epoch (a std::string) and an
 * EpochBus<std::string> as the epoch source.
 *----------------------------------------------------------------------------*/
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <poll.h>
#include <sys/socket.h>
#include <unistd.h>
#include <gtest/gtest.h>
#include "nrtk/ntrip_front.hpp"

using namespace nrtk;

namespace {

using steady=std::chrono::steady_clock;

double since(steady::time_point t)
{
    return std::chrono::duration<double>(steady::now()-t).count();
}
bool wait_for(const std::function<bool()> &pred, double timeout_s=3.0)
{
    auto t0=steady::now();
    while (!pred()) {
        if (since(t0)>timeout_s) return false;
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    return true;
}
/* "$" + body + "*hh" */
std::string nmea(const std::string &body)
{
    unsigned cs=0;
    for (char c:body) cs^=(unsigned char)c;
    char t[8];
    snprintf(t,sizeof(t),"*%02X",cs);
    return "$"+body+t;
}
bool starts(const std::string &s, const std::string &prefix)
{
    return s.compare(0,prefix.size(),prefix)==0;
}
std::string body_of(const std::string &resp)
{
    size_t e=resp.find("\r\n\r\n");
    return e==std::string::npos?std::string():resp.substr(e+4);
}
long content_length(const std::string &resp)
{
    size_t p=resp.find("\r\nContent-Length: ");
    return p==std::string::npos?-1:atol(resp.c_str()+p+18);
}
const char GGA_BODY[]   ="GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,46.9,M,,";
const char GGA_Q0_BODY[]="GPGGA,123519,4807.038,N,01131.000,E,0,08,0.9,545.4,M,46.9,M,,";

/* raw TCP client ---------------------------------------------------------------*/
class Client {
public:
    explicit Client(int port) {
        sockaddr_in a;
        fd_=socket(AF_INET,SOCK_STREAM,0);
#ifdef SO_NOSIGPIPE
        int on=1;
        setsockopt(fd_,SOL_SOCKET,SO_NOSIGPIPE,&on,sizeof(on));
#endif
        memset(&a,0,sizeof(a));
        a.sin_family=AF_INET;
        a.sin_port=htons((uint16_t)port);
        a.sin_addr.s_addr=htonl(INADDR_LOOPBACK);
        if (connect(fd_,(sockaddr *)&a,sizeof(a))<0) {close(fd_); fd_=-1;}
    }
    ~Client() {shut();}
    Client(const Client &)=delete;
    Client &operator=(const Client &)=delete;
    bool ok() const {return fd_>=0;}
    void shut() {if (fd_>=0) close(fd_); fd_=-1;}
    bool send(const std::string &s) {
#ifdef MSG_NOSIGNAL
        const int flags=MSG_NOSIGNAL;
#else
        const int flags=0;
#endif
        for (size_t off=0;off<s.size();) {
            ssize_t n=::send(fd_,s.data()+off,s.size()-off,flags);
            if (n<=0) return false;
            off+=(size_t)n;
        }
        return true;
    }
    /* receive until done() (true), EOF or timeout (false) */
    bool fill(double timeout_s, const std::function<bool()> &done) {
        auto t0=steady::now();
        while (!done()) {
            if (eof||since(t0)>timeout_s) return false;
            pollfd p={fd_,POLLIN,0};
            if (poll(&p,1,10)<=0) continue;
            char b[4096];
            ssize_t n=recv(fd_,b,sizeof(b),0);
            if (n<=0) eof=true;
            else buf.append(b,(size_t)n);
        }
        return true;
    }
    bool read_until(const std::string &delim, double timeout_s=3.0) {
        return fill(timeout_s,[&]{return buf.find(delim)!=std::string::npos;});
    }
    bool read_n(size_t n, double timeout_s=3.0) {
        return fill(timeout_s,[&]{return buf.size()>=n;});
    }
    std::string take(size_t n) {
        std::string s=buf.substr(0,n);
        buf.erase(0,n);
        return s;
    }
    /* everything up to EOF (true: EOF seen) */
    bool read_all(double timeout_s=3.0) {
        fill(timeout_s,[&]{return eof;});
        return eof;
    }
    /* response header block */
    bool read_head(std::string &head, double timeout_s=3.0) {
        if (!read_until("\r\n\r\n",timeout_s)) return false;
        head=take(buf.find("\r\n\r\n")+4);
        return true;
    }
    /* one HTTP chunk (empty data: last chunk) */
    bool read_chunk(std::string &data, double timeout_s=3.0) {
        if (!read_until("\r\n",timeout_s)) return false;
        size_t e=buf.find("\r\n");
        size_t n=strtoul(buf.substr(0,e).c_str(),nullptr,16);
        buf.erase(0,e+2);
        if (!read_n(n+2,timeout_s)||buf.compare(n,2,"\r\n")) return false;
        data=take(n);
        buf.erase(0,2);
        return true;
    }
    std::string buf;
    bool eof=false;
private:
    int fd_=-1;
};

/* fake session: echoes the epoch string, "CLOSE" ends the stream -----------------*/
struct Shared {
    std::mutex mtx;
    std::vector<GgaFix> gga;
    std::vector<ClientInfo> info;
    std::atomic<int> created{0},destroyed{0};
    std::atomic<bool> refuse{false},overlap{false};
    size_t ngga() {std::lock_guard<std::mutex> lk(mtx); return gga.size();}
};
class FakeSession : public StreamSession {
public:
    FakeSession(std::shared_ptr<Shared> sh, const ClientInfo &ci) : sh_(std::move(sh)), ci_(ci) {
        sh_->created++;
    }
    ~FakeSession() override {sh_->destroyed++;}
    void on_gga(const GgaFix &fix) override {
        Busy b(this);
        std::lock_guard<std::mutex> lk(sh_->mtx);
        sh_->gga.push_back(fix);
    }
    bool on_epoch(const void *epoch, std::vector<uint8_t> &out) override {
        Busy b(this);
        const std::string &s=*static_cast<const std::string *>(epoch);
        out.assign(s.begin(),s.end());
        return s!="CLOSE";
    }
    std::string status_json() const override {
        Busy b(this);
        return "{\"fake\":"+std::to_string(ci_.id)+"}";
    }
private:
    struct Busy {           /* flags concurrent calls on one session */
        const FakeSession *s;
        explicit Busy(const FakeSession *p) : s(p) {
            if (s->busy_.fetch_add(1)) s->sh_->overlap=true;
            std::this_thread::sleep_for(std::chrono::microseconds(200));
        }
        ~Busy() {s->busy_.fetch_sub(1);}
    };
    std::shared_ptr<Shared> sh_;
    ClientInfo ci_;
    mutable std::atomic<int> busy_{0};
};

} // namespace

/* GGA parser -----------------------------------------------------------------------*/
TEST(GgaParse, ValidSentence)
{
    std::string s=nmea(GGA_BODY);
    EXPECT_EQ(s.substr(s.size()-3),"*47");      /* well-known example checksum */
    GgaFix f;
    ASSERT_TRUE(parse_gga(s+"\r\n",f));
    EXPECT_TRUE(f.valid);
    EXPECT_EQ(f.quality,1);
    EXPECT_EQ(f.nsat,8);
    EXPECT_DOUBLE_EQ(f.hdop,0.9);
    EXPECT_NEAR(f.tod,12*3600+35*60+19,1E-9);
    double lat=48.0+7.038/60.0,lon=11.0+31.0/60.0;
    EXPECT_NEAR(f.lat*R2D,lat,1E-12);
    EXPECT_NEAR(f.lon*R2D,lon,1E-12);
    EXPECT_DOUBLE_EQ(f.h_msl,545.4);
    EXPECT_TRUE(f.geoid_from_gga);
    EXPECT_DOUBLE_EQ(f.geoid,46.9);
    EXPECT_NEAR(f.h_ell,592.3,1E-9);
    double pos[3]={lat*D2R,lon*D2R,592.3},r[3],llh[3];
    pos2ecef(pos,r);
    for (int i=0;i<3;i++) EXPECT_NEAR(f.ecef[i],r[i],1E-6);
    ecef2pos(f.ecef,llh);
    EXPECT_NEAR(llh[0]*R2D,lat,1E-9);
    EXPECT_NEAR(llh[1]*R2D,lon,1E-9);
    EXPECT_NEAR(llh[2],592.3,1E-4);

    /* without checksum; lower-case checksum, blanks around, fractional time */
    ASSERT_TRUE(parse_gga(std::string("$")+GGA_BODY,f));
    EXPECT_TRUE(f.valid);
    ASSERT_TRUE(parse_gga(" $GPGGA,000001.5,0000.0000,N,00000.0000,E,1,04,1.0,0.0,M,0.0,M,,*6d \n",f));
    EXPECT_TRUE(f.valid);
    EXPECT_NEAR(f.tod,1.5,1E-12);
}

TEST(GgaParse, ChecksumMismatchRejected)
{
    GgaFix f;
    std::string s=nmea(GGA_BODY);
    ASSERT_TRUE(parse_gga(s,f));
    s[s.size()-1]='8';                          /* *48 */
    EXPECT_FALSE(parse_gga(s,f));
    EXPECT_FALSE(f.valid);
    EXPECT_FALSE(parse_gga(s.substr(0,s.size()-1),f));         /* "*4" */
    std::string t=nmea(GGA_BODY);
    t[10]='6';                                  /* corrupted payload */
    EXPECT_FALSE(parse_gga(t,f));
}

TEST(GgaParse, QualityZeroNotValid)
{
    GgaFix f;
    ASSERT_TRUE(parse_gga(nmea(GGA_Q0_BODY),f));
    EXPECT_FALSE(f.valid);
    EXPECT_EQ(f.quality,0);
    ASSERT_TRUE(parse_gga("$GPGGA,,,,,,0,00,99.99,,,,,,*48",f));      /* receiver without fix */
    EXPECT_FALSE(f.valid);
    EXPECT_EQ(f.ecef[0],0.0);
    ASSERT_TRUE(parse_gga(nmea("GPGGA,123519,,,,,1,08,0.9,545.4,M,46.9,M,,"),f));  /* no position */
    EXPECT_FALSE(f.valid);
}

TEST(GgaParse, MissingGeoidUsesEgm96)
{
    GgaFix f;
    ASSERT_TRUE(parse_gga(nmea("GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,,M,,"),f));
    EXPECT_TRUE(f.valid);
    EXPECT_FALSE(f.geoid_from_gga);
    double pos[3]={f.lat,f.lon,0.0};
    double n=geoidh(pos);
    EXPECT_GT(n,40.0);                          /* EGM96 near Munich: ~47 m */
    EXPECT_LT(n,55.0);
    EXPECT_DOUBLE_EQ(f.geoid,n);
    EXPECT_DOUBLE_EQ(f.h_ell,545.4+n);
    double r[3];
    pos[2]=545.4+n;
    pos2ecef(pos,r);
    for (int i=0;i<3;i++) EXPECT_NEAR(f.ecef[i],r[i],1E-6);

    /* sentence truncated after the altitude */
    ASSERT_TRUE(parse_gga(nmea("GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4"),f));
    EXPECT_FALSE(f.geoid_from_gga);
    EXPECT_DOUBLE_EQ(f.geoid,n);
}

TEST(GgaParse, MultiConstellationTalkers)
{
    GgaFix f;
    ASSERT_TRUE(parse_gga(nmea("GNGGA,001043.00,3347.5123,S,15112.3456,W,4,12,0.6,10.0,M,-20.1,M,1.0,0000"),f));
    EXPECT_TRUE(f.valid);
    EXPECT_EQ(f.quality,4);
    EXPECT_EQ(f.nsat,12);
    EXPECT_NEAR(f.tod,643.0,1E-9);
    EXPECT_NEAR(f.lat*R2D,-(33.0+47.5123/60.0),1E-12);
    EXPECT_NEAR(f.lon*R2D,-(151.0+12.3456/60.0),1E-12);
    EXPECT_NEAR(f.h_ell,-10.1,1E-9);
    for (const char *t:{"GL","GA","GB","BD","GQ"}) {
        std::string b=std::string(t)+(GGA_BODY+2);
        EXPECT_TRUE(parse_gga(nmea(b),f)) << b;
        EXPECT_TRUE(f.valid) << b;
    }
}

TEST(GgaParse, RejectsMalformed)
{
    GgaFix f;
    EXPECT_FALSE(parse_gga("",f));
    EXPECT_FALSE(parse_gga(GGA_BODY,f));                                        /* no '$' */
    EXPECT_FALSE(parse_gga(nmea("GPRMC,123519,A,4807.038,N,01131.000,E,022.4,084.4,230394,003.1,W"),f));
    EXPECT_FALSE(parse_gga(nmea("GPGSA,A,3,04,05,,09,12,,,24,,,,,2.5,1.3,2.1"),f));
    EXPECT_FALSE(parse_gga(nmea("GPGGA,123519,4807.038,N,01131.000"),f));        /* too short */
    EXPECT_FALSE(parse_gga(nmea("GPGGA,123519,48x7.038,N,01131.000,E,1,08,0.9,545.4,M,46.9,M,,"),f));
    EXPECT_FALSE(parse_gga(nmea("GPGGA,123519,4867.038,N,01131.000,E,1,08,0.9,545.4,M,46.9,M,,"),f));
    EXPECT_FALSE(parse_gga(nmea("GPGGA,123519,4807.038,X,01131.000,E,1,08,0.9,545.4,M,46.9,M,,"),f));
    EXPECT_FALSE(parse_gga(nmea("GPGGA,123519,4807.038,N,18131.000,E,1,08,0.9,545.4,M,46.9,M,,"),f));
    EXPECT_FALSE(parse_gga(nmea("GPGGA,253519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,46.9,M,,"),f));
    EXPECT_FALSE(parse_gga(nmea("GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,F,46.9,M,,"),f));
    EXPECT_FALSE(parse_gga(nmea("GPGGA,123519,4807.038,N,01131.000,E,q,08,0.9,545.4,M,46.9,M,,"),f));
    EXPECT_FALSE(parse_gga(nmea("GPGGA,123519,4807.038,N,01131.000,E,1,08,0.9,5e2,M,46.9,M,,"),f));
    EXPECT_FALSE(parse_gga(nmea(GGA_BODY)+"xx",f));
}

/* base64 --------------------------------------------------------------------------*/
TEST(NtripBase64, RoundTrip)
{
    const char *v[][2]={{"",""},{"f","Zg=="},{"fo","Zm8="},{"foo","Zm9v"},{"foob","Zm9vYg=="},
                        {"fooba","Zm9vYmE="},{"foobar","Zm9vYmFy"},{"user:pass","dXNlcjpwYXNz"}};
    for (auto &t:v) {
        EXPECT_EQ(base64_encode(t[0]),t[1]);
        EXPECT_EQ(base64_decode(t[1]),t[0]);
    }
    std::string bin;
    for (int i=0;i<300;i++) {
        bin+=(char)(i*37&0xFF);
        ASSERT_EQ(base64_decode(base64_encode(bin)),bin) << i;
    }
    EXPECT_EQ(base64_decode("Zm9v\r\nYmE"),"fooba");     /* whitespace, no padding */
    EXPECT_EQ(base64_decode("Zm9v*mE="),"");             /* invalid character */
}

/* caster ---------------------------------------------------------------------------*/
class NtripServerTest : public ::testing::Test {
protected:
    void SetUp() override {
        opt.port=0;
        opt.bind_addr="127.0.0.1";
        opt.caster_host="caster.test";
        MountPoint m;
        m.name="VRS"; m.identifier="Taipei"; m.format="RTCM 3.3";
        m.format_details="1006(10),1077(1)"; m.nav_system="GPS+GAL"; m.network="NRTK";
        m.country="TWN"; m.lat=25.03; m.lon=121.56; m.need_gga=true;
        mounts.push_back(m);
        m.name="NOGGA"; m.need_gga=false;
        mounts.push_back(m);
    }
    void TearDown() override {
        if (srv) srv->stop();
        EXPECT_FALSE(sh->overlap.load()) << "session methods called concurrently";
        EXPECT_EQ(sh->created.load(),sh->destroyed.load());
    }
    bool start(HttpHandler http=nullptr) {
        std::shared_ptr<Shared> s=sh;
        srv.reset(new NtripServer(opt,mounts,
            [s](const ClientInfo &ci)->std::unique_ptr<StreamSession> {
                {std::lock_guard<std::mutex> lk(s->mtx); s->info.push_back(ci);}
                if (s->refuse) return nullptr;
                return std::unique_ptr<StreamSession>(new FakeSession(s,ci));
            },
            [this](uint64_t *seq, int timeout_ms)->std::shared_ptr<const void> {
                return bus.wait(seq,timeout_ms);
            },
            std::move(http)));
        if (!srv->start()) return false;
        port=srv->port();
        return port>0;
    }
    void publish(const std::string &s) {bus.publish(std::make_shared<std::string>(s));}
    /* NTRIP v1 mountpoint request: hdr = extra header lines, data = bytes sent
     * right after the header block. true on "ICY 200 OK" */
    bool open_v1(Client &c, const std::string &mount, const std::string &hdr="",
                 const std::string &data="") {
        if (!c.ok()||!c.send("GET /"+mount+" HTTP/1.0\r\nUser-Agent: NTRIP test/1.0\r\n"+hdr+"\r\n"+data)) {
            return false;
        }
        if (!c.fill(3.0,[&]{return c.buf.size()>=12||c.buf.find("\r\n\r\n")!=std::string::npos;})) {
            return false;
        }
        if (starts(c.buf,"ICY 200 OK\r\n")) {c.take(12); return true;}
        c.read_all();
        head=c.buf;
        return false;
    }
    /* NTRIP v2 mountpoint request. true on "HTTP/1.1 200 OK", header in head */
    bool open_v2(Client &c, const std::string &mount, const std::string &hdr="") {
        head.clear();
        return c.ok()&&c.send("GET /"+mount+" HTTP/1.1\r\nHost: localhost\r\nNtrip-Version: Ntrip/2.0\r\n"
                              "User-Agent: NTRIP test/2.0\r\n"+hdr+"\r\n")&&
               c.read_head(head)&&starts(head,"HTTP/1.1 200 OK\r\n");
    }
    NtripServerOptions opt;
    std::vector<MountPoint> mounts;
    std::shared_ptr<Shared> sh=std::make_shared<Shared>();
    EpochBus<std::string> bus;
    std::unique_ptr<NtripServer> srv;
    int port=0;
    std::string head;
};

TEST_F(NtripServerTest, SourcetableV1)
{
    ASSERT_TRUE(start());
    Client c(port);
    ASSERT_TRUE(c.send("GET / HTTP/1.0\r\nUser-Agent: NTRIP test/1.0\r\n\r\n"));
    ASSERT_TRUE(c.read_all());
    const std::string &r=c.buf;
    EXPECT_TRUE(starts(r,"SOURCETABLE 200 OK\r\nServer: NTRIP nrtk/")) << r;
    EXPECT_NE(r.find("\r\nContent-Type: text/plain\r\n"),std::string::npos);
    std::string body=body_of(r);
    EXPECT_EQ(content_length(r),(long)body.size());
    EXPECT_EQ(body,srv->sourcetable());
    EXPECT_NE(body.find("STR;VRS;Taipei;RTCM 3.3;1006(10),1077(1);2;GPS+GAL;NRTK;TWN;25.03;121.56;"
                        "1;1;nrtk;none;N;N;9600;none\r\n"),std::string::npos) << body;
    EXPECT_NE(body.find("STR;NOGGA;Taipei;RTCM 3.3;1006(10),1077(1);2;GPS+GAL;NRTK;TWN;25.03;121.56;"
                        "0;1;nrtk;none;N;N;9600;none\r\n"),std::string::npos) << body;
    EXPECT_NE(body.find("CAS;caster.test;"+std::to_string(port)+";NRTK;nrtk;0;TWN;25.03;121.56;"),
              std::string::npos) << body;
    ASSERT_GE(body.size(),16u);
    EXPECT_EQ(body.substr(body.size()-16),"ENDSOURCETABLE\r\n");
}

TEST_F(NtripServerTest, SourcetableV2)
{
    ASSERT_TRUE(start());
    Client c(port);
    ASSERT_TRUE(c.send("GET / HTTP/1.1\nHost: localhost\nNtrip-Version: Ntrip/2.0\n"
                       "User-Agent: NTRIP test/2.0\n\n"));                  /* LF only */
    ASSERT_TRUE(c.read_all());
    const std::string &r=c.buf;
    EXPECT_TRUE(starts(r,"HTTP/1.1 200 OK\r\n")) << r;
    EXPECT_NE(r.find("\r\nNtrip-Version: Ntrip/2.0\r\n"),std::string::npos);
    EXPECT_NE(r.find("\r\nContent-Type: gnss/sourcetable\r\n"),std::string::npos);
    std::string body=body_of(r);
    EXPECT_EQ(content_length(r),(long)body.size());
    EXPECT_NE(body.find("STR;VRS;Taipei;"),std::string::npos);
    ASSERT_GE(body.size(),16u);
    EXPECT_EQ(body.substr(body.size()-16),"ENDSOURCETABLE\r\n");
}

TEST_F(NtripServerTest, UnknownPathNotFound)
{
    ASSERT_TRUE(start());
    Client c(port);
    ASSERT_TRUE(c.send("GET /NOSUCH HTTP/1.0\r\nUser-Agent: NTRIP test/1.0\r\n\r\n"));
    ASSERT_TRUE(c.read_all());
    EXPECT_TRUE(starts(c.buf,"HTTP/1.1 404 Not Found\r\n")) << c.buf;
    EXPECT_EQ(sh->created.load(),0);
}

TEST_F(NtripServerTest, HttpHandler)
{
    std::vector<std::string> paths;
    std::mutex pm;
    ASSERT_TRUE(start([&](const std::string &path, std::string &ctype, std::string &body) {
        {std::lock_guard<std::mutex> lk(pm); paths.push_back(path);}
        if (path=="/boom") throw std::runtime_error("boom");
        if (path!="/status") return false;
        ctype="application/json";
        body="{\"clients\":"+std::to_string(srv->clients().size())+"}";  /* re-entrant */
        return true;
    }));
    Client s(port);
    ASSERT_TRUE(open_v1(s,"NOGGA")) << head;
    {
        Client c(port);
        ASSERT_TRUE(c.send("GET /status?pretty=1 HTTP/1.1\r\nHost: localhost\r\n\r\n"));
        ASSERT_TRUE(c.read_all());
        EXPECT_TRUE(starts(c.buf,"HTTP/1.1 200 OK\r\n")) << c.buf;
        EXPECT_NE(c.buf.find("\r\nContent-Type: application/json\r\n"),std::string::npos);
        EXPECT_EQ(body_of(c.buf),"{\"clients\":1}");
        EXPECT_EQ(content_length(c.buf),13);
    }
    {
        Client c(port);
        ASSERT_TRUE(c.send("GET /metrics2 HTTP/1.1\r\nHost: localhost\r\n\r\n"));
        ASSERT_TRUE(c.read_all());
        EXPECT_TRUE(starts(c.buf,"HTTP/1.1 404 Not Found\r\n")) << c.buf;
    }
    {
        Client c(port);
        ASSERT_TRUE(c.send("GET /boom HTTP/1.1\r\nHost: localhost\r\n\r\n"));
        ASSERT_TRUE(c.read_all());
        EXPECT_TRUE(starts(c.buf,"HTTP/1.1 500 Internal Server Error\r\n")) << c.buf;
    }
    std::lock_guard<std::mutex> lk(pm);
    ASSERT_EQ(paths.size(),3u);
    EXPECT_EQ(paths[0],"/status");
    EXPECT_EQ(paths[1],"/metrics2");
}

TEST_F(NtripServerTest, BadRequests)
{
    ASSERT_TRUE(start());
    {
        Client c(port);
        ASSERT_TRUE(c.send("POST /VRS HTTP/1.1\r\nHost: localhost\r\nNtrip-Version: Ntrip/2.0\r\n\r\n"));
        ASSERT_TRUE(c.read_all());
        EXPECT_TRUE(starts(c.buf,"HTTP/1.1 405 Method Not Allowed\r\n")) << c.buf;
    }
    {
        Client c(port);
        ASSERT_TRUE(c.send("HELLO\r\n\r\n"));
        ASSERT_TRUE(c.read_all());
        EXPECT_TRUE(starts(c.buf,"HTTP/1.1 400 Bad Request\r\n")) << c.buf;
    }
    {
        Client c(port);                         /* header block > 8 KB */
        std::string h="GET /VRS HTTP/1.0\r\n";
        while (h.size()<9000) h+="X-Pad: 0123456789012345678901234567890123456789\r\n";
        ASSERT_TRUE(c.send(h));
        ASSERT_TRUE(c.read_all());
        EXPECT_TRUE(starts(c.buf,"HTTP/1.1 400 Bad Request\r\n")) << c.buf;
    }
    EXPECT_EQ(sh->created.load(),0);
}

TEST_F(NtripServerTest, GgaInsteadOfBlankLine)
{
    ASSERT_TRUE(start());
    Client a(port),b(port);
    /* GGA before the blank line */
    ASSERT_TRUE(open_v1(a,"VRS",nmea(GGA_BODY)+"\r\n")) << head;
    /* GGA and no blank line at all (control characters in the agent are replaced) */
    ASSERT_TRUE(b.send("GET /VRS HTTP/1.0\r\nUser-Agent: NTRIP\x01x\tv\r\n"+nmea(GGA_BODY)+"\r\n"));
    ASSERT_TRUE(b.read_n(12));
    EXPECT_EQ(b.take(12),"ICY 200 OK\r\n");
    ASSERT_TRUE(wait_for([&]{return sh->ngga()==2;}));
    publish("X");
    ASSERT_TRUE(a.read_n(1));
    ASSERT_TRUE(b.read_n(1));
    std::lock_guard<std::mutex> lk(sh->mtx);
    ASSERT_EQ(sh->info.size(),2u);
    EXPECT_EQ(sh->info[1].agent,"NTRIP?x v");
}

TEST_F(NtripServerTest, StreamV1)
{
    ASSERT_TRUE(start());
    Client c(port);
    /* GGA in the same segment as the request header */
    ASSERT_TRUE(open_v1(c,"VRS","",nmea(GGA_BODY)+"\r\n")) << head;
    ASSERT_TRUE(wait_for([&]{return sh->ngga()==1;}));
    publish("EPOCH-1");
    ASSERT_TRUE(c.read_n(7));
    EXPECT_EQ(c.take(7),"EPOCH-1");

    /* garbage, other NMEA, a GGA split over two writes, quality 0 (ignored) */
    std::string g=nmea("GNGGA,123520,4807.038,N,01131.000,E,4,12,0.6,545.4,M,46.9,M,1.0,0000");
    ASSERT_TRUE(c.send("\x01\x02garbage\r\n"+nmea("GPGSA,A,3,04,05,,09,12,,,24,,,,,2.5,1.3,2.1")+
                       "\r\n"+g.substr(0,20)));
    std::this_thread::sleep_for(std::chrono::milliseconds(60));
    ASSERT_TRUE(c.send(g.substr(20)+"\r\n"+nmea(GGA_Q0_BODY)+"\r\n"));
    ASSERT_TRUE(wait_for([&]{return sh->ngga()==2;}));
    publish("EPOCH-2");
    ASSERT_TRUE(c.read_n(7));
    EXPECT_EQ(c.take(7),"EPOCH-2");
    {
        std::lock_guard<std::mutex> lk(sh->mtx);
        ASSERT_EQ(sh->gga.size(),2u);
        EXPECT_EQ(sh->gga[1].quality,4);
        EXPECT_NEAR(sh->gga[0].lat*R2D,48.0+7.038/60.0,1E-12);
        ASSERT_EQ(sh->info.size(),1u);
        EXPECT_EQ(sh->info[0].mount,"VRS");
        EXPECT_EQ(sh->info[0].ntrip_version,1);
        EXPECT_EQ(sh->info[0].agent,"NTRIP test/1.0");
        EXPECT_EQ(sh->info[0].user,"");
        EXPECT_TRUE(starts(sh->info[0].peer,"127.0.0.1:")) << sh->info[0].peer;
    }
    auto cl=srv->clients();
    ASSERT_EQ(cl.size(),1u);
    EXPECT_EQ(cl[0].info.mount,"VRS");
    EXPECT_EQ(cl[0].ngga,2);
    EXPECT_EQ(cl[0].bytes_sent,14);
    EXPECT_GT(cl[0].connected_s,0.0);
    EXPECT_EQ(cl[0].session_json,"{\"fake\":"+std::to_string(cl[0].info.id)+"}");

    /* status polling concurrent with epochs and GGA */
    std::atomic<bool> polling{true};
    std::thread th([&]{
        while (polling) {
            srv->clients();
            std::this_thread::sleep_for(std::chrono::microseconds(100));
        }
    });
    for (int i=0;i<5;i++) {
        std::string e="E"+std::to_string(i);
        ASSERT_TRUE(c.send(nmea(GGA_BODY)+"\r\n"));
        publish(e);
        ASSERT_TRUE(c.read_n(e.size()));
        EXPECT_EQ(c.take(e.size()),e);
    }
    polling=false;
    th.join();

    /* session ends the stream */
    publish("CLOSE");
    ASSERT_TRUE(c.read_all());
    EXPECT_EQ(c.buf,"CLOSE");
    EXPECT_TRUE(wait_for([&]{return sh->destroyed.load()==1&&srv->clients().empty();}));
}

TEST_F(NtripServerTest, StreamV2Chunked)
{
    ASSERT_TRUE(start());
    Client c(port);
    ASSERT_TRUE(open_v2(c,"VRS","Ntrip-GGA: "+nmea(GGA_BODY)+"\r\n")) << head;
    EXPECT_NE(head.find("\r\nNtrip-Version: Ntrip/2.0\r\n"),std::string::npos);
    EXPECT_NE(head.find("\r\nContent-Type: gnss/data\r\n"),std::string::npos);
    EXPECT_NE(head.find("\r\nTransfer-Encoding: chunked\r\n"),std::string::npos);
    EXPECT_NE(head.find("\r\nCache-Control: no-store, no-cache, max-age=0\r\n"),std::string::npos);
    ASSERT_TRUE(wait_for([&]{return sh->ngga()==1;}));      /* from the Ntrip-GGA header */

    std::string big(1000,'x'),d;
    for (const std::string &e:{std::string("ABC"),std::string("HELLO WORLD"),big}) {
        publish(e);
        ASSERT_TRUE(c.read_chunk(d));
        EXPECT_EQ(d,e);
    }
    ASSERT_TRUE(c.send(nmea(GGA_BODY)+"\r\n"));             /* upstream GGA */
    ASSERT_TRUE(wait_for([&]{return sh->ngga()==2;}));
    auto cl=srv->clients();
    ASSERT_EQ(cl.size(),1u);
    EXPECT_EQ(cl[0].info.ntrip_version,2);
    EXPECT_EQ(cl[0].info.agent,"NTRIP test/2.0");
    EXPECT_EQ(cl[0].bytes_sent,(long)(3+11+big.size()));

    publish("CLOSE");
    ASSERT_TRUE(c.read_chunk(d));
    EXPECT_EQ(d,"CLOSE");
    ASSERT_TRUE(c.read_chunk(d));                           /* last chunk */
    EXPECT_EQ(d,"");
    EXPECT_TRUE(c.read_all());
    EXPECT_TRUE(c.buf.empty());
}

TEST_F(NtripServerTest, BasicAuth)
{
    UserAccount u;
    u.user="alice"; u.pass="s3cret"; u.max_conn=1;
    opt.users.push_back(u);
    opt.allow_anonymous=false;
    ASSERT_TRUE(start());
    EXPECT_NE(srv->sourcetable().find(";25.03;121.56;1;1;nrtk;none;B;N;9600;none\r\n"),
              std::string::npos);
    std::string good="Authorization: Basic "+base64_encode("alice:s3cret")+"\r\n";
    std::string bad="Authorization: Basic "+base64_encode("alice:wrong")+"\r\n";
    {
        Client c(port);                         /* v1, no credentials */
        EXPECT_FALSE(open_v1(c,"VRS"));
        EXPECT_TRUE(starts(head,"HTTP/1.0 401 Unauthorized\r\n")) << head;
    }
    {
        Client c(port);                         /* v1, wrong password */
        EXPECT_FALSE(open_v1(c,"VRS",bad));
        EXPECT_TRUE(starts(head,"HTTP/1.0 401 Unauthorized\r\n")) << head;
    }
    {
        Client c(port);                         /* v2, wrong password */
        EXPECT_FALSE(open_v2(c,"VRS",bad));
        EXPECT_TRUE(starts(head,"HTTP/1.1 401 Unauthorized\r\n")) << head;
        EXPECT_NE(head.find("\r\nWWW-Authenticate: Basic realm=\"/VRS\"\r\n"),std::string::npos);
    }
    EXPECT_EQ(sh->created.load(),0);

    Client a(port);                             /* correct password */
    ASSERT_TRUE(open_v1(a,"VRS",good,nmea(GGA_BODY)+"\r\n")) << head;
    auto cl=srv->clients();
    ASSERT_EQ(cl.size(),1u);
    EXPECT_EQ(cl[0].info.user,"alice");
    {
        Client b(port);                         /* account max_conn=1 */
        EXPECT_FALSE(open_v2(b,"VRS",good));
        EXPECT_TRUE(starts(head,"HTTP/1.1 503 Service Unavailable\r\n")) << head;
    }
    a.shut();                                   /* slot released on close */
    bool ok=false;
    for (int i=0;i<40&&!ok;i++) {
        Client t(port);
        ok=open_v1(t,"NOGGA",good);
        if (!ok) std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    EXPECT_TRUE(ok);
}

TEST_F(NtripServerTest, GgaTimeout)
{
    opt.gga_timeout_s=0.5;
    ASSERT_TRUE(start());
    Client a(port),b(port),c(port),d(port);
    ASSERT_TRUE(open_v1(a,"VRS"));                          /* no GGA */
    ASSERT_TRUE(open_v2(b,"VRS","Ntrip-GGA: "+nmea(GGA_Q0_BODY)+"\r\n"));  /* invalid GGA */
    ASSERT_TRUE(open_v1(c,"NOGGA"));                        /* mount without GGA */
    ASSERT_TRUE(open_v1(d,"VRS","",nmea(GGA_BODY)+"\r\n")); /* valid GGA */
    auto t0=steady::now();
    ASSERT_TRUE(a.read_all(3.0));
    ASSERT_TRUE(b.read_all(3.0));
    double dt=since(t0);
    EXPECT_GT(dt,0.3);
    EXPECT_LT(dt,2.0);
    EXPECT_EQ(b.buf,"0\r\n\r\n");                           /* v2 last chunk */
    std::this_thread::sleep_for(std::chrono::milliseconds(300));
    EXPECT_FALSE(c.read_all(0.05));                         /* still open */
    EXPECT_FALSE(d.read_all(0.05));
    publish("E1");
    ASSERT_TRUE(c.read_n(2));
    ASSERT_TRUE(d.read_n(2));
    EXPECT_EQ(c.take(2),"E1");
    EXPECT_EQ(d.take(2),"E1");
    EXPECT_EQ(srv->clients().size(),2u);
}

TEST_F(NtripServerTest, MaxClients)
{
    opt.max_clients=1;
    ASSERT_TRUE(start());
    Client a(port);
    ASSERT_TRUE(open_v1(a,"NOGGA"));
    {
        Client b(port);
        EXPECT_FALSE(open_v1(b,"NOGGA"));
        EXPECT_TRUE(starts(head,"HTTP/1.0 503 Service Unavailable\r\n")) << head;
    }
    {
        Client s(port);                         /* sourcetable is not a stream */
        ASSERT_TRUE(s.send("GET / HTTP/1.0\r\n\r\n"));
        ASSERT_TRUE(s.read_all());
        EXPECT_TRUE(starts(s.buf,"SOURCETABLE 200 OK\r\n"));
    }
    a.shut();
    bool ok=false;
    for (int i=0;i<40&&!ok;i++) {
        Client t(port);
        ok=open_v1(t,"NOGGA");
        if (!ok) std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    EXPECT_TRUE(ok);
}

TEST_F(NtripServerTest, FactoryRefusal)
{
    sh->refuse=true;
    ASSERT_TRUE(start());
    Client c(port);
    EXPECT_FALSE(open_v2(c,"VRS"));
    EXPECT_TRUE(starts(head,"HTTP/1.1 503 Service Unavailable\r\n")) << head;
    {
        std::lock_guard<std::mutex> lk(sh->mtx);
        EXPECT_EQ(sh->info.size(),1u);
    }
    EXPECT_TRUE(srv->clients().empty());
}

TEST_F(NtripServerTest, StopWithClients)
{
    ASSERT_TRUE(start());
    Client a(port),b(port),c(port),d(port);
    ASSERT_TRUE(open_v1(a,"VRS","",nmea(GGA_BODY)+"\r\n"));
    ASSERT_TRUE(open_v2(b,"NOGGA"));
    ASSERT_TRUE(c.send("GET /VRS HTTP/1.0\r\n"));           /* incomplete header */
    ASSERT_TRUE(d.ok());                                     /* connected, silent */
    ASSERT_TRUE(wait_for([&]{return srv->clients().size()==2;}));
    auto t0=steady::now();
    srv->stop();
    double dt=since(t0);
    EXPECT_LT(dt,1.0);
    EXPECT_TRUE(srv->clients().empty());
    EXPECT_EQ(sh->destroyed.load(),2);
    EXPECT_TRUE(a.read_all(1.0));
    EXPECT_TRUE(b.read_all(1.0));
    EXPECT_EQ(b.buf,"0\r\n\r\n");
    EXPECT_TRUE(c.read_all(1.0));
    EXPECT_TRUE(d.read_all(1.0));
    srv->stop();                                             /* idempotent */
    Client e(port);                                          /* listener closed */
    EXPECT_FALSE(e.ok());
}
