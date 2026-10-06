/*------------------------------------------------------------------------------
 * ntrip_server.cpp : NTRIP caster front end (v1 and v2) on POSIX sockets
 *
 * threads: one accept thread + one thread per TCP connection. A connection
 * thread reads the request header, then either answers (sourcetable, HTTP
 * handler, error) and closes, or runs the stream loop of one StreamSession:
 *     poll(socket, 50 ms) -> upstream NMEA -> parse_gga -> session->on_gga()
 *     wait_epoch(&seq, 100 ms) -> session->on_epoch() -> send (v2: chunked)
 * Session methods are called only by the connection thread, or by clients()
 * (status_json), always under the per-connection mutex: never concurrently.
 * Finished connection threads are joined by the accept thread; stop() joins
 * the accept thread and every connection thread.
 *
 * replies: 400 bad request, 401 missing/wrong credentials (anonymous access
 * disabled), 404 unknown path, 405 method other than GET, 500 HTTP handler
 * exception, 503 max_clients reached, account max_conn reached or the session
 * factory returned null.
 *----------------------------------------------------------------------------*/
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <list>
#include <map>
#include <string>
#include <vector>
#include <fcntl.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <poll.h>
#include <pthread.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>
#include "nrtk/common.h"
#include "nrtk/ntrip_front.hpp"

namespace nrtk {

using steady=std::chrono::steady_clock;

namespace {

const size_t HDR_MAX       =8192;       /* request header block limit (bytes) */
const double HDR_TIMEOUT_S =5.0;        /* request header read timeout (s) */
const double SEND_TIMEOUT_S=5.0;        /* close if the client takes no data (s) */
const double LINGER_S      =1.0;        /* drain input before close (s) */
const int    POLL_MS       =50;         /* socket poll interval (ms) */
const int    EPOCH_WAIT_MS =100;        /* wait_epoch timeout (ms) */
const int    ACCEPT_POLL_MS=100;        /* accept poll interval (ms) */
const size_t NMEA_MAX      =1024;       /* upstream sentence length limit */
const int    EXTRA_CONN    =64;         /* connections beyond max_clients (HTTP, handshakes) */
const size_t THREAD_STACK  =8u<<20;     /* connection thread stack (macOS default: 512 KB) */
const char   SERVER_ID[]   ="NTRIP nrtk/" NRTK_VERSION;

#ifdef MSG_NOSIGNAL
const int SEND_FLAGS=MSG_NOSIGNAL;
#else
const int SEND_FLAGS=0;
#endif

double secs(steady::duration d)
{
    return std::chrono::duration<double>(d).count();
}
std::string lower(std::string s)
{
    for (auto &c:s) c=(char)tolower((unsigned char)c);
    return s;
}
std::string trim(const std::string &s)
{
    size_t a=s.find_first_not_of(" \t"),b=s.find_last_not_of(" \t");
    return a==std::string::npos?std::string():s.substr(a,b-a+1);
}
/* sourcetable field: no separators or line breaks */
std::string field(const std::string &s)
{
    std::string o=s;
    for (auto &c:o) if (c==';') c=','; else if (c=='\r'||c=='\n') c=' ';
    return o;
}
std::string fmt2(double v)
{
    char s[32];
    snprintf(s,sizeof(s),"%.2f",v);
    return s;
}
/* compare credentials without early exit */
bool same(const std::string &a, const std::string &b)
{
    if (a.size()!=b.size()) return false;
    unsigned char d=0;
    for (size_t i=0;i<a.size();i++) d|=(unsigned char)(a[i]^b[i]);
    return d==0;
}

/* sockets ------------------------------------------------------------------------*/
void set_nonblock(int fd)
{
    fcntl(fd,F_SETFL,fcntl(fd,F_GETFL,0)|O_NONBLOCK);
    fcntl(fd,F_SETFD,FD_CLOEXEC);
}
void setup_client(int fd)
{
    int on=1;
    set_nonblock(fd);
#ifdef SO_NOSIGPIPE
    setsockopt(fd,SOL_SOCKET,SO_NOSIGPIPE,&on,sizeof(on));
#endif
    setsockopt(fd,IPPROTO_TCP,TCP_NODELAY,&on,sizeof(on));
}
std::string peer_str(const sockaddr_storage &sa, socklen_t len)
{
    char host[256]="?",serv[32]="0";
    getnameinfo((const sockaddr *)&sa,len,host,sizeof(host),serv,sizeof(serv),
                NI_NUMERICHOST|NI_NUMERICSERV);
    std::string h=host;
    if (h.compare(0,7,"::ffff:")==0&&h.find('.')!=std::string::npos) h.erase(0,7);
    if (h.find(':')!=std::string::npos) h="["+h+"]";
    return h+":"+serv;
}
/* send all bytes on a non-blocking socket. false: error, stop requested or the
 * client took no data for SEND_TIMEOUT_S */
bool send_all(int fd, const void *data, size_t n, const std::atomic<bool> &run)
{
    const char *p=(const char *)data;
    auto last=steady::now();

    while (n>0) {
        ssize_t k=send(fd,p,n,SEND_FLAGS);
        if (k>0) {p+=k; n-=(size_t)k; last=steady::now(); continue;}
        if (k<0&&errno==EINTR) continue;
        if (k<0&&(errno==EAGAIN||errno==EWOULDBLOCK)) {
            if (!run||secs(steady::now()-last)>SEND_TIMEOUT_S) return false;
            pollfd pf={fd,POLLOUT,0};
            poll(&pf,1,POLL_MS);
            continue;
        }
        return false;
    }
    return true;
}
/* close: optionally half-close and drain input first, so that a client still
 * sending (e.g. GGA) does not get a reset that destroys our last reply */
void close_socket(int fd, bool linger, const std::atomic<bool> &run)
{
    if (fd<0) return;
    if (linger&&run) {
        char buf[1024];
        auto t0=steady::now();
        shutdown(fd,SHUT_WR);
        while (run&&secs(steady::now()-t0)<LINGER_S) {
            pollfd pf={fd,POLLIN,0};
            int r=poll(&pf,1,POLL_MS);
            if (r<0&&errno!=EINTR) break;
            if (r<=0) continue;
            ssize_t k=recv(fd,buf,sizeof(buf),0);
            if (k==0||(k<0&&errno!=EAGAIN&&errno!=EWOULDBLOCK&&errno!=EINTR)) break;
        }
    }
    close(fd);
}

/* threads with an explicit stack size (session code may need a deep stack) ------*/
void *thread_entry(void *arg)
{
    std::unique_ptr<std::function<void()>> f((std::function<void()> *)arg);
    try {(*f)();}
    catch (...) {}
    return nullptr;
}
bool start_thread(pthread_t *th, std::function<void()> fn)
{
    pthread_attr_t at;
    auto *f=new std::function<void()>(std::move(fn));

    pthread_attr_init(&at);
    pthread_attr_setstacksize(&at,THREAD_STACK);
    int r=pthread_create(th,&at,thread_entry,f);
    pthread_attr_destroy(&at);
    if (r) {delete f; return false;}
    return true;
}

/* HTTP request ------------------------------------------------------------------*/
struct Request {
    std::string method,target;
    std::map<std::string,std::string> hdr;      /* lower-case names */
    std::string get(const char *name) const {
        auto it=hdr.find(name);
        return it==hdr.end()?std::string():it->second;
    }
};
/* end of the header block: after an empty line (CRLF or LF), or at a line
 * starting with '$' (rovers sending GGA without the blank line). npos: more */
size_t header_end(const std::string &b)
{
    size_t p=b.find_first_not_of("\r\n");
    if (p==std::string::npos) return std::string::npos;
    for (bool first=true;;first=false) {
        if (!first&&p<b.size()&&b[p]=='$') return p;
        size_t nl=b.find('\n',p);
        if (nl==std::string::npos) return std::string::npos;
        if (nl==p||(nl==p+1&&b[p]=='\r')) return nl+1;
        p=nl+1;
    }
}
/* read the header block. 1: ok, 0: timeout/closed/stop, -1: too large.
 * bytes received after the header block are returned in rest */
int read_header(int fd, const std::atomic<bool> &run, std::string &hdr, std::string &rest)
{
    std::string b;
    char buf[2048];
    auto t0=steady::now();

    for (;;) {
        size_t e=header_end(b);
        if (e!=std::string::npos) {hdr=b.substr(0,e); rest=b.substr(e); return 1;}
        if (b.size()>HDR_MAX) return -1;
        if (!run||secs(steady::now()-t0)>HDR_TIMEOUT_S) return 0;
        pollfd pf={fd,POLLIN,0};
        int r=poll(&pf,1,POLL_MS);
        if (r<0&&errno!=EINTR) return 0;
        if (r<=0) continue;
        ssize_t k=recv(fd,buf,sizeof(buf),0);
        if (k==0) return 0;
        if (k<0) {
            if (errno==EAGAIN||errno==EWOULDBLOCK||errno==EINTR) continue;
            return 0;
        }
        b.append(buf,(size_t)k);
    }
}
bool parse_request(const std::string &h, Request &q)
{
    bool first=true;
    for (size_t p=0;p<h.size();) {
        size_t nl=h.find('\n',p);
        std::string line=h.substr(p,nl==std::string::npos?std::string::npos:nl-p);
        p=nl==std::string::npos?h.size():nl+1;
        if (!line.empty()&&line.back()=='\r') line.pop_back();
        if (line.empty()) {if (first) continue; break;}
        if (line[0]=='$') break;
        for (auto &ch:line) {           /* no control characters in names/values (logs, JSON) */
            if (ch=='\t') ch=' ';
            else if ((unsigned char)ch<0x20||ch==0x7f) ch='?';
        }
        if (first) {                            /* METHOD SP TARGET [SP VERSION] */
            std::vector<std::string> tok;
            for (size_t i=0;i<line.size();) {
                size_t j=line.find(' ',i);
                if (j==std::string::npos) j=line.size();
                if (j>i) tok.push_back(line.substr(i,j-i));
                i=j+1;
            }
            if (tok.size()<2) return false;
            q.method=tok[0];
            q.target=tok.size()==2&&tok[1].compare(0,5,"HTTP/")==0?"/":tok[1];
            first=false;
            continue;
        }
        size_t c=line.find(':');
        if (c==std::string::npos) continue;
        q.hdr[lower(trim(line.substr(0,c)))]=trim(line.substr(c+1));
    }
    return !first;
}
/* path of the request target: absolute URI reduced to its path, no query */
std::string req_path(const std::string &target)
{
    std::string p=target;
    if (lower(p.substr(0,7))=="http://") {
        size_t s=p.find('/',7);
        p=s==std::string::npos?"/":p.substr(s);
    }
    size_t q=p.find('?');
    if (q!=std::string::npos) p.erase(q);
    return p;
}

/* replies --------------------------------------------------------------------------*/
const char *reason_phrase(int code)
{
    switch (code) {
        case 200: return "OK";
        case 400: return "Bad Request";
        case 401: return "Unauthorized";
        case 404: return "Not Found";
        case 405: return "Method Not Allowed";
        case 503: return "Service Unavailable";
    }
    return "Internal Server Error";
}
/* HTTP reply (always closes the connection) */
std::string http_reply(bool http10, bool v2, int code, const std::string &ctype,
                       const std::string &body, const std::string &extra=std::string())
{
    std::string s=std::string(http10?"HTTP/1.0 ":"HTTP/1.1 ")+std::to_string(code)+" "+
                  reason_phrase(code)+"\r\n";
    if (v2) s+="Ntrip-Version: Ntrip/2.0\r\n";
    s+=std::string("Server: ")+SERVER_ID+"\r\n"+extra;
    s+="Content-Type: "+ctype+"\r\nContent-Length: "+std::to_string(body.size())+
       "\r\nConnection: close\r\n\r\n"+body;
    return s;
}

/* upstream NMEA: split received bytes into sentences. A sentence starts at '$'
 * and ends at CR/LF, at the next '$' or right after a "*hh" checksum. Complete
 * sentences go to out, the incomplete tail stays in buf. Other bytes are
 * ignored. */
void nmea_split(std::string &buf, const char *p, size_t n, std::vector<std::string> &out)
{
    for (size_t i=0;i<n;i++) {
        char c=p[i];
        if (c=='\r'||c=='\n'||c=='$') {
            if (!buf.empty()) out.push_back(buf);
            buf.clear();
            if (c=='$') buf+=c;
            continue;
        }
        if (buf.empty()) continue;
        if (buf.size()>=NMEA_MAX) {buf.clear(); continue;}
        buf+=c;
        size_t m=buf.size();
        if (m>=4&&buf[m-3]=='*'&&isxdigit((unsigned char)buf[m-2])&&
            isxdigit((unsigned char)buf[m-1])) {
            out.push_back(buf);
            buf.clear();
        }
    }
}

} // namespace

/* implementation -----------------------------------------------------------------*/
struct NtripServer::Impl {
    struct Conn {
        uint64_t id=0;
        int fd=-1;
        std::string peer;
        pthread_t th{};
        bool th_ok=false;               /* thread started (accept thread / stop only) */
        std::atomic<bool> done{false};  /* thread finished, ready to join */
        steady::time_point t0;
        mutable std::mutex mtx;         /* session calls and the fields below */
        ClientInfo info;
        std::unique_ptr<StreamSession> session;
        long bytes=0;
        int ngga=0;
    };
    /* stream slot (max_clients, account max_conn), released on destruction */
    struct Slot {
        Impl *m=nullptr;
        std::string user;               /* account name, empty: anonymous */
        ~Slot() {release();}
        void release() {
            if (!m) return;
            std::lock_guard<std::mutex> lk(m->mtx);
            m->nstream--;
            if (!user.empty()&&--m->nuser[user]<=0) m->nuser.erase(user);
            m=nullptr;
        }
    };

    NtripServerOptions opt;
    std::vector<MountPoint> mounts;
    SessionFactory factory;
    std::function<std::shared_ptr<const void>(uint64_t *,int)> wait_epoch;
    HttpHandler http;

    std::mutex ctl;                     /* start/stop */
    std::atomic<bool> run{false};
    std::atomic<int> bport{0};
    int lfd=-1;
    std::thread acc_th;                 /* accept thread */

    mutable std::mutex mtx;             /* conns and counters */
    std::list<std::shared_ptr<Conn>> conns;
    uint64_t next_id=1;
    int nstream=0;
    std::map<std::string,int> nuser;

    std::string table() const;
    void accept_loop();
    void reap();
    void conn_main(const std::shared_ptr<Conn> &c);
    void serve(Conn &c);
    void reply(Conn &c, const std::string &resp, int code, const std::string &mount,
               const char *why);
    void stream(Conn &c, const MountPoint &mp, bool v2, const std::string &gga_hdr,
                const std::string &rest, Slot &slot);
};

/* sourcetable --------------------------------------------------------------------*/
std::string NtripServer::Impl::table() const
{
    bool auth=!opt.users.empty()&&!opt.allow_anonymous;
    int port=bport>0?bport.load():opt.port;
    std::string s;
    const MountPoint *m0=mounts.empty()?nullptr:&mounts[0];

    s+="CAS;"+field(opt.caster_host)+";"+std::to_string(port)+";"+field(opt.network_name)+";"+
       field(opt.operator_name)+";0;"+(m0?field(m0->country):std::string())+";"+
       fmt2(m0?m0->lat:0.0)+";"+fmt2(m0?m0->lon:0.0)+";0.0.0.0;0;none\r\n";
    for (const MountPoint &m:mounts) {
        s+="STR;"+field(m.name)+";"+field(m.identifier)+";"+field(m.format)+";"+
           field(m.format_details)+";2;"+field(m.nav_system)+";"+
           field(m.network.empty()?opt.network_name:m.network)+";"+field(m.country)+";"+
           fmt2(m.lat)+";"+fmt2(m.lon)+";"+(m.need_gga?"1":"0")+";1;"+field(m.generator)+
           ";none;"+(auth?"B":"N")+";N;9600;none\r\n";
    }
    s+="ENDSOURCETABLE\r\n";
    return s;
}

/* accept thread --------------------------------------------------------------------*/
void NtripServer::Impl::accept_loop()
{
    int maxconn=(opt.max_clients>0?opt.max_clients:1024)+EXTRA_CONN;
    steady::time_point last_warn;

    while (run) {
        reap();
        pollfd pf={lfd,POLLIN,0};
        int r=poll(&pf,1,ACCEPT_POLL_MS);
        if (r<0&&errno!=EINTR) {
            std::this_thread::sleep_for(std::chrono::milliseconds(ACCEPT_POLL_MS));
            continue;
        }
        if (r<=0||!run) continue;
        sockaddr_storage sa;
        socklen_t sl=sizeof(sa);
        int fd=accept(lfd,(sockaddr *)&sa,&sl);
        if (fd<0) {                     /* e.g. EMFILE: back off, warn every 10 s */
            if (errno!=EAGAIN&&errno!=EWOULDBLOCK&&errno!=EINTR&&errno!=ECONNABORTED) {
                if (secs(steady::now()-last_warn)>10.0) {
                    nrtk_log(NRTK_LOG_WARN,"ntrip","accept: %s",strerror(errno));
                    last_warn=steady::now();
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(ACCEPT_POLL_MS));
            }
            continue;
        }
        setup_client(fd);
        auto c=std::make_shared<Conn>();
        c->fd=fd;
        c->peer=peer_str(sa,sl);
        c->t0=steady::now();
        bool full=false;
        {
            std::lock_guard<std::mutex> lk(mtx);
            c->id=next_id++;
            if ((int)conns.size()>=maxconn) full=true;
            else conns.push_back(c);
        }
        if (full) {
            nrtk_log(NRTK_LOG_WARN,"ntrip","reject id=%llu peer=%s: too many connections (%d)",
                     (unsigned long long)c->id,c->peer.c_str(),maxconn);
            close(fd);
            continue;
        }
        nrtk_log(NRTK_LOG_DEBUG,"ntrip","accept id=%llu peer=%s",(unsigned long long)c->id,
                 c->peer.c_str());
        if (start_thread(&c->th,[this,c]{conn_main(c);})) {c->th_ok=true; continue;}
        nrtk_log(NRTK_LOG_ERROR,"ntrip","id=%llu: cannot start connection thread",
                 (unsigned long long)c->id);
        {std::lock_guard<std::mutex> lk(mtx); conns.remove(c);}
        close(fd);
    }
}
/* join finished connection threads */
void NtripServer::Impl::reap()
{
    std::vector<std::shared_ptr<Conn>> dead;
    {
        std::lock_guard<std::mutex> lk(mtx);
        for (auto it=conns.begin();it!=conns.end();) {
            if ((*it)->done) {dead.push_back(*it); it=conns.erase(it);}
            else ++it;
        }
    }
    for (auto &c:dead) if (c->th_ok) pthread_join(c->th,nullptr);
}

/* connection thread ------------------------------------------------------------------*/
void NtripServer::Impl::conn_main(const std::shared_ptr<Conn> &c)
{
    try {
        serve(*c);
    }
    catch (const std::exception &e) {
        nrtk_log(NRTK_LOG_ERROR,"ntrip","id=%llu exception: %s",(unsigned long long)c->id,e.what());
    }
    catch (...) {
        nrtk_log(NRTK_LOG_ERROR,"ntrip","id=%llu unknown exception",(unsigned long long)c->id);
    }
    std::unique_ptr<StreamSession> s;          /* after an exception in stream() */
    {std::lock_guard<std::mutex> lk(c->mtx); s=std::move(c->session);}
    s.reset();
    if (c->fd>=0) {close(c->fd); c->fd=-1;}
    c->done=true;
}
/* send a reply (HTTP or error) and close */
void NtripServer::Impl::reply(Conn &c, const std::string &resp, int code,
                              const std::string &mount, const char *why)
{
    if (code>=400) {
        nrtk_log(NRTK_LOG_INFO,"ntrip","reject id=%llu peer=%s mount=%s: %d %s",
                 (unsigned long long)c.id,c.peer.c_str(),mount.c_str(),code,why);
    }
    else {
        nrtk_log(NRTK_LOG_DEBUG,"ntrip","http id=%llu peer=%s %s: %d (%zu bytes)",
                 (unsigned long long)c.id,c.peer.c_str(),mount.c_str(),code,resp.size());
    }
    bool ok=send_all(c.fd,resp.data(),resp.size(),run);
    close_socket(c.fd,ok,run);
    c.fd=-1;
}
void NtripServer::Impl::serve(Conn &c)
{
    std::string hdr,rest;
    Request q;

    int r=read_header(c.fd,run,hdr,rest);
    if (r<0) {
        reply(c,http_reply(false,false,400,"text/plain","request header too large\r\n"),400,"",
              "header too large");
        return;
    }
    if (r==0) {
        nrtk_log(NRTK_LOG_DEBUG,"ntrip","close id=%llu peer=%s: no request (timeout or closed)",
                 (unsigned long long)c.id,c.peer.c_str());
        close_socket(c.fd,false,run);
        c.fd=-1;
        return;
    }
    if (!parse_request(hdr,q)) {
        reply(c,http_reply(false,false,400,"text/plain","bad request\r\n"),400,"","bad request line");
        return;
    }
    bool v2=lower(q.get("ntrip-version")).compare(0,7,"ntrip/2")==0;
    std::string path=req_path(q.target);
    std::string mount=path.size()>0&&path[0]=='/'?path.substr(1):path;

    if (q.method!="GET") {
        reply(c,http_reply(false,v2,405,"text/plain","method not allowed\r\n","Allow: GET\r\n"),
              405,mount,("method "+q.method).c_str());
        return;
    }
    /* sourcetable */
    if (mount.empty()) {
        std::string body=table(),resp;
        if (v2) resp=http_reply(false,true,200,"gnss/sourcetable",body);
        else {
            resp=std::string("SOURCETABLE 200 OK\r\nServer: ")+SERVER_ID+
                 "\r\nContent-Type: text/plain\r\nContent-Length: "+std::to_string(body.size())+
                 "\r\n\r\n"+body;
        }
        reply(c,resp,200,"/","sourcetable");
        return;
    }
    const MountPoint *mp=nullptr;
    for (auto &m:mounts) if (m.name==mount) {mp=&m; break;}

    /* other paths: HTTP handler (/status, /metrics ...) */
    if (!mp) {
        std::string ctype="text/plain",body;
        bool ok=false;
        if (http) {
            std::string err;
            try {
                ok=http(path,ctype,body);
            }
            catch (const std::exception &e) {err=e.what();}
            catch (...) {err="unknown exception";}
            if (!err.empty()) {
                nrtk_log(NRTK_LOG_ERROR,"ntrip","http handler %s: %s",path.c_str(),err.c_str());
                reply(c,http_reply(false,v2,500,"text/plain","internal error\r\n"),500,mount,
                      "handler exception");
                return;
            }
        }
        if (ok) reply(c,http_reply(false,false,200,ctype,body),200,path,"http");
        else reply(c,http_reply(false,v2,404,"text/plain","not found\r\n"),404,mount,"not found");
        return;
    }
    /* authentication (HTTP basic) */
    std::string user,pass,auth=q.get("authorization");
    bool cred=false;
    if (lower(auth.substr(0,6))=="basic ") {
        std::string d=base64_decode(trim(auth.substr(6)));
        size_t k=d.find(':');
        if (k!=std::string::npos) {user=d.substr(0,k); pass=d.substr(k+1); cred=true;}
    }
    const UserAccount *acc=nullptr;
    if (cred) for (auto &u:opt.users) if (u.user==user&&same(u.pass,pass)) {acc=&u; break;}
    if (!acc&&!opt.allow_anonymous) {
        reply(c,http_reply(!v2,v2,401,"text/plain","401 Unauthorized\r\n",
                           "WWW-Authenticate: Basic realm=\"/"+mount+"\"\r\n"),
              401,mount,cred?"bad credentials":"no credentials");
        return;
    }
    /* stream slot */
    Slot slot;
    const char *full=nullptr;
    {
        std::lock_guard<std::mutex> lk(mtx);
        auto it=acc?nuser.find(acc->user):nuser.end();
        if (opt.max_clients>0&&nstream>=opt.max_clients) full="max_clients reached";
        else if (acc&&acc->max_conn>0&&it!=nuser.end()&&it->second>=acc->max_conn) {
            full="account max_conn reached";
        }
        else {
            nstream++;
            if (acc) nuser[acc->user]++;
            slot.m=this;
            slot.user=acc?acc->user:"";
        }
    }
    if (full) {
        reply(c,http_reply(!v2,v2,503,"text/plain","503 Service Unavailable\r\n"),503,mount,full);
        return;
    }
    /* session */
    ClientInfo ci;
    ci.id=c.id;
    ci.mount=mp->name;
    ci.user=acc?acc->user:(opt.users.empty()?user:"");
    ci.peer=c.peer;
    ci.agent=q.get("user-agent");
    ci.ntrip_version=v2?2:1;
    std::unique_ptr<StreamSession> s;
    if (factory) {
        std::string err;
        try {
            s=factory(ci);
        }
        catch (const std::exception &e) {err=e.what();}
        catch (...) {err="unknown exception";}
        if (!err.empty()) {
            s.reset();
            nrtk_log(NRTK_LOG_ERROR,"ntrip","id=%llu session factory: %s",(unsigned long long)c.id,
                     err.c_str());
        }
    }
    if (!s) {
        slot.release();
        reply(c,http_reply(!v2,v2,503,"text/plain","503 Service Unavailable\r\n"),503,mount,
              "no session");
        return;
    }
    {
        std::lock_guard<std::mutex> lk(c.mtx);
        c.info=ci;
        c.session=std::move(s);
    }
    nrtk_log(NRTK_LOG_INFO,"ntrip","connect id=%llu peer=%s mount=%s user=%s agent=\"%s\" ntrip=v%d",
             (unsigned long long)c.id,c.peer.c_str(),ci.mount.c_str(),ci.user.c_str(),
             ci.agent.c_str(),ci.ntrip_version);
    stream(c,*mp,v2,q.get("ntrip-gga"),rest,slot);
}

/* stream loop ------------------------------------------------------------------------*/
void NtripServer::Impl::stream(Conn &c, const MountPoint &mp, bool v2, const std::string &gga_hdr,
                               const std::string &rest, Slot &slot)
{
    std::string nbuf;
    std::vector<std::string> lines;
    bool have_gga=false,linger=true,fin=true;
    const char *why=nullptr;
    char buf[4096];
    uint64_t seq=0;

    /* response header */
    std::string h=v2?std::string("HTTP/1.1 200 OK\r\nNtrip-Version: Ntrip/2.0\r\nServer: ")+SERVER_ID+
                     "\r\nContent-Type: gnss/data\r\nCache-Control: no-store, no-cache, max-age=0\r\n"
                     "Transfer-Encoding: chunked\r\nConnection: close\r\n\r\n"
                    :std::string("ICY 200 OK\r\n");   /* as RTKLIB and most v1 rovers expect */
    if (!send_all(c.fd,h.data(),h.size(),run)) {why="send failed"; linger=fin=false;}
    auto t_start=steady::now();

    /* upstream bytes -> NMEA sentences -> valid GGA to the session */
    auto upstream=[&](const char *p, size_t n){
        lines.clear();
        nmea_split(nbuf,p,n,lines);
        for (auto &l:lines) {
            GgaFix f;
            if (!parse_gga(l,f)||!f.valid) continue;
            {
                std::lock_guard<std::mutex> lk(c.mtx);
                c.session->on_gga(f);
                c.ngga++;
            }
            if (!have_gga) {
                nrtk_log(NRTK_LOG_INFO,"ntrip","gga id=%llu lat=%.6f lon=%.6f h=%.1f q=%d nsat=%d",
                         (unsigned long long)c.id,f.lat*R2D,f.lon*R2D,f.h_ell,f.quality,f.nsat);
            }
            have_gga=true;
        }
    };
    if (!why) {
        if (!gga_hdr.empty()) upstream((gga_hdr+"\r\n").c_str(),gga_hdr.size()+2);
        if (!rest.empty()) upstream(rest.data(),rest.size());
    }
    while (!why) {
        if (!run) {why="server stop"; linger=false; break;}
        if (mp.need_gga&&!have_gga&&opt.gga_timeout_s>0.0&&
            secs(steady::now()-t_start)>opt.gga_timeout_s) {
            why="no valid GGA";
            break;
        }
        pollfd pf={c.fd,POLLIN,0};
        int r=poll(&pf,1,POLL_MS);
        if (r<0&&errno!=EINTR) {why="poll error"; linger=fin=false; break;}
        if (r>0) {
            ssize_t k=recv(c.fd,buf,sizeof(buf),0);
            if (k==0) {why="client closed"; linger=fin=false; break;}
            if (k<0&&errno!=EAGAIN&&errno!=EWOULDBLOCK&&errno!=EINTR) {
                why="recv error"; linger=fin=false; break;
            }
            if (k>0) upstream(buf,(size_t)k);
        }
        if (!wait_epoch) continue;
        std::shared_ptr<const void> ep=wait_epoch(&seq,EPOCH_WAIT_MS);
        if (!ep) continue;
        std::vector<uint8_t> out;
        bool keep;
        {
            std::lock_guard<std::mutex> lk(c.mtx);
            keep=c.session->on_epoch(ep.get(),out);
        }
        if (!out.empty()) {
            bool ok;
            if (v2) {                           /* HTTP chunk */
                char ch[24];
                int nh=snprintf(ch,sizeof(ch),"%zx\r\n",out.size());
                std::string b(ch,(size_t)nh);
                b.append((const char *)out.data(),out.size());
                b+="\r\n";
                ok=send_all(c.fd,b.data(),b.size(),run);
            }
            else ok=send_all(c.fd,out.data(),out.size(),run);
            if (!ok) {
                if (!run) {why="server stop"; linger=false;}
                else {why="send failed or timeout"; linger=fin=false;}
                break;
            }
            std::lock_guard<std::mutex> lk(c.mtx);
            c.bytes+=(long)out.size();
        }
        if (!keep) why="session end";
    }
    /* end: session (no more calls), slot (by the caller), then the socket */
    std::unique_ptr<StreamSession> s;
    long bytes;
    int ngga;
    {
        std::lock_guard<std::mutex> lk(c.mtx);
        s=std::move(c.session);
        bytes=c.bytes;
        ngga=c.ngga;
    }
    nrtk_log(NRTK_LOG_INFO,"ntrip","close id=%llu peer=%s mount=%s: %s, %ld bytes, %d GGA, %.1f s",
             (unsigned long long)c.id,c.peer.c_str(),mp.name.c_str(),why,bytes,ngga,
             secs(steady::now()-c.t0));
    s.reset();
    slot.release();
    if (fin&&v2) (void)send(c.fd,"0\r\n\r\n",5,SEND_FLAGS);    /* last chunk, best effort */
    close_socket(c.fd,linger,run);
    c.fd=-1;
}

/* NtripServer ------------------------------------------------------------------------*/
NtripServer::NtripServer(const NtripServerOptions &opt, std::vector<MountPoint> mounts,
                         SessionFactory factory,
                         std::function<std::shared_ptr<const void>(uint64_t *seq,int timeout_ms)> wait_epoch,
                         HttpHandler http)
    : impl_(new Impl)
{
    impl_->opt=opt;
    impl_->mounts=std::move(mounts);
    impl_->factory=std::move(factory);
    impl_->wait_epoch=std::move(wait_epoch);
    impl_->http=std::move(http);
}
NtripServer::~NtripServer()
{
    stop();
}
bool NtripServer::start()
{
    Impl &m=*impl_;
    std::lock_guard<std::mutex> lk(m.ctl);
    if (m.run) return true;
    if (m.opt.port<0||m.opt.port>65535) {
        nrtk_log(NRTK_LOG_ERROR,"ntrip","invalid port %d",m.opt.port);
        return false;
    }
    std::string host=m.opt.bind_addr.empty()?"0.0.0.0":m.opt.bind_addr;
    std::string svc=std::to_string(m.opt.port);
    addrinfo hints,*res=nullptr;
    memset(&hints,0,sizeof(hints));
    hints.ai_family=AF_UNSPEC;
    hints.ai_socktype=SOCK_STREAM;
    hints.ai_flags=AI_PASSIVE;
    int err=getaddrinfo(host.c_str(),svc.c_str(),&hints,&res);
    if (err) {
        nrtk_log(NRTK_LOG_ERROR,"ntrip","bind address %s: %s",host.c_str(),gai_strerror(err));
        return false;
    }
    int fd=-1,lasterr=0,on=1;
    for (addrinfo *ai=res;ai;ai=ai->ai_next) {
        if ((fd=socket(ai->ai_family,ai->ai_socktype,ai->ai_protocol))<0) {lasterr=errno; continue;}
        setsockopt(fd,SOL_SOCKET,SO_REUSEADDR,&on,sizeof(on));
        if (bind(fd,ai->ai_addr,ai->ai_addrlen)==0&&listen(fd,SOMAXCONN)==0) break;
        lasterr=errno;
        close(fd);
        fd=-1;
    }
    freeaddrinfo(res);
    if (fd<0) {
        nrtk_log(NRTK_LOG_ERROR,"ntrip","cannot listen on %s:%d: %s",host.c_str(),m.opt.port,
                 strerror(lasterr));
        return false;
    }
    set_nonblock(fd);
    sockaddr_storage sa;
    socklen_t sl=sizeof(sa);
    int port=m.opt.port;
    if (getsockname(fd,(sockaddr *)&sa,&sl)==0) {
        if (sa.ss_family==AF_INET) port=ntohs(((sockaddr_in *)&sa)->sin_port);
        else if (sa.ss_family==AF_INET6) port=ntohs(((sockaddr_in6 *)&sa)->sin6_port);
    }
    m.lfd=fd;
    m.bport=port;
    m.run=true;
    try {
        m.acc_th=std::thread(&Impl::accept_loop,&m);
    }
    catch (const std::exception &e) {
        nrtk_log(NRTK_LOG_ERROR,"ntrip","cannot start accept thread: %s",e.what());
        m.run=false;
        close(m.lfd);
        m.lfd=-1;
        return false;
    }
    if (m.opt.users.empty()&&!m.opt.allow_anonymous) {
        nrtk_log(NRTK_LOG_WARN,"ntrip","no accounts and anonymous access disabled: "
                 "all stream requests will be rejected");
    }
    nrtk_log(NRTK_LOG_INFO,"ntrip","listening on %s:%d, %zu mountpoints, %zu accounts%s",
             host.c_str(),port,m.mounts.size(),m.opt.users.size(),
             m.opt.allow_anonymous?", anonymous allowed":"");
    return true;
}
void NtripServer::stop()
{
    Impl &m=*impl_;
    std::lock_guard<std::mutex> lk(m.ctl);
    if (!m.run&&!m.acc_th.joinable()) return;
    m.run=false;
    if (m.acc_th.joinable()) m.acc_th.join();
    if (m.lfd>=0) {close(m.lfd); m.lfd=-1;}
    std::list<std::shared_ptr<Impl::Conn>> v;
    {
        std::lock_guard<std::mutex> lk2(m.mtx);
        v.swap(m.conns);
    }
    for (auto &c:v) if (c->th_ok) pthread_join(c->th,nullptr);
    nrtk_log(NRTK_LOG_INFO,"ntrip","stopped (%zu connections closed)",v.size());
}
int NtripServer::port() const
{
    return impl_->bport>0?impl_->bport.load():impl_->opt.port;
}
std::vector<NtripServer::ClientStatus> NtripServer::clients() const
{
    std::vector<std::shared_ptr<Impl::Conn>> v;
    std::vector<ClientStatus> out;
    {
        std::lock_guard<std::mutex> lk(impl_->mtx);
        v.assign(impl_->conns.begin(),impl_->conns.end());
    }
    auto now=steady::now();
    for (auto &c:v) {
        std::lock_guard<std::mutex> lk(c->mtx);
        if (!c->session) continue;
        ClientStatus s;
        s.info=c->info;
        s.connected_s=secs(now-c->t0);
        s.bytes_sent=c->bytes;
        s.ngga=c->ngga;
        try {
            s.session_json=c->session->status_json();
        }
        catch (...) {}
        out.push_back(std::move(s));
    }
    return out;
}
std::string NtripServer::sourcetable() const
{
    return impl_->table();
}

} // namespace nrtk
