/*------------------------------------------------------------------------------
 * offline.cpp : offline input (RINEX) and VRS output writers
 *----------------------------------------------------------------------------*/
#include <cstring>
#include <sys/stat.h>
#include "nrtk/pipeline.hpp"
#include "nrtk/rtcm_out.h"

namespace nrtk {

int mkdirs(const std::string &dir)
{
    for (size_t i=1;i<=dir.size();i++) {
        if (i==dir.size()||dir[i]=='/') mkdir(dir.substr(0,i).c_str(),0755);
    }
    struct stat st;
    return stat(dir.c_str(),&st)==0&&S_ISDIR(st.st_mode);
}

/* offline reader ----------------------------------------------------------------*/
OfflineReader::OfflineReader(const nrtk_config_t &cfg) : cfg_(cfg)
{
    src_.resize(cfg.nsta);
    open_.assign(cfg.nsta,0);
    memset(&nav_,0,sizeof(nav_));
    int nopen=0;
    for (int i=0;i<cfg.nsta;i++) {
        if (!cfg.sta[i].enabled||!cfg.sta[i].obsfile[0]) continue;
        if (nrtk_rnxsrc_open(&src_[i],cfg.sta[i].obsfile,i+1)) {open_[i]=1; nopen++;}
    }
    if (!nrtk_nav_read(cfg.nav,&nav_)) {
        nrtk_log(NRTK_LOG_ERROR,"offline","no navigation data (offline.nav)");
        return;
    }
    tint_=cfg.tint>0.0?cfg.tint:1.0;
    t_=cfg.ts; te_=cfg.te;
    if (!t_.time) {   /* earliest first epoch, rounded to the grid */
        for (int i=0;i<cfg.nsta;i++) {
            if (!open_[i]) continue;
            gtime_t t=nrtk_rnxsrc_peek(&src_[i]);
            if (t.time&&(!t_.time||timediff(t,t_)<0.0)) t_=t;
        }
        double tow=time2gpst(t_,nullptr);
        t_=timeadd(t_,ceil(tow/tint_-1E-6)*tint_-tow);
    }
    ok_=nopen>=2&&t_.time;
    nrtk_log(NRTK_LOG_INFO,"offline","%d observation files, start %s",nopen,time_str(t_).c_str());
}
OfflineReader::~OfflineReader()
{
    for (int i=0;i<cfg_.nsta;i++) if (open_[i]) nrtk_rnxsrc_close(&src_[i]);
    nrtk_nav_free(&nav_);
}
bool OfflineReader::next(gtime_t &t, std::vector<std::vector<obsd_t>> &raw)
{
    if (!ok_) return false;
    if (!first_) t_=timeadd(t_,tint_);
    first_=false;
    if (te_.time&&timediff(t_,te_)>DTTOL) return false;
    raw.assign(cfg_.nsta,std::vector<obsd_t>());
    int alive=0;
    std::vector<obsd_t> buf(MAXOBS);
    for (int i=0;i<cfg_.nsta;i++) {
        if (!open_[i]) continue;
        int n=nrtk_rnxsrc_read(&src_[i],t_,buf.data(),MAXOBS);
        if (n<0) continue;
        alive++;
        if (n>0) raw[i].assign(buf.begin(),buf.begin()+n);
    }
    t=t_;
    return alive>0;
}

/* VRS writers -------------------------------------------------------------------*/
struct VrsWriter::RtcmState {
    nrtk_rtcmout_t enc;
    int mode=0, sta=-1;
};

VrsWriter::VrsWriter(const nrtk_config_t &cfg, const std::string &dir,
                     const std::string &name, const double *pos, int staid,
                     const nav_t *nav) : cfg_(cfg)
{
    memcpy(pos_,pos,sizeof(pos_));
    mkdirs(dir);
    std::string base=dir+"/"+name;
    rnx_=fopen((base+".obs").c_str(),"w");
    rtcm_=fopen((base+".rtcm3").c_str(),"wb");
    if (cfg.diag_csv) {
        diag_=fopen((base+"_diag.csv").c_str(),"w");
        epo_=fopen((base+"_epoch.csv").c_str(),"w");
        if (diag_) fprintf(diag_,"time,sat,master,isref,ok,nbl,dd_iono,dd_nondisp,dd_trop_model,"
                                 "rms_iono,rms_nd,a_iono,b_iono,a_nd,b_nd,extrap,newfix\n");
        if (epo_) fprintf(epo_,"time,mode,master,nsat,nok,nlockreset,iono_activity,trop_resid,"
                                "degraded,proc_ms\n");
    }
    memset(&ropt_,0,sizeof(ropt_));
    ropt_.rnxver=304;
    ropt_.navsys=SYS_GPS|SYS_GAL|SYS_CMP|SYS_QZS;
    ropt_.obstype=OBSTYPE_PR|OBSTYPE_CP|OBSTYPE_DOP|OBSTYPE_SNR;
    ropt_.freqtype=FREQTYPE_ALL;
    ropt_.tint=cfg.tint;
    for (int m=0;m<RNX_NUMSYS;m++) memset(ropt_.mask[m],'1',MAXCODE);
    snprintf(ropt_.prog,sizeof(ropt_.prog),"nrtk-gen %s",NRTK_VERSION);
    snprintf(ropt_.runby,sizeof(ropt_.runby),"nrtk");
    snprintf(ropt_.marker,sizeof(ropt_.marker),"%s",name.c_str());
    snprintf(ropt_.markertype,sizeof(ropt_.markertype),"NON_PHYSICAL");
    snprintf(ropt_.rec[1],sizeof(ropt_.rec[1]),"NRTK NETWORK RTK");
    snprintf(ropt_.rec[2],sizeof(ropt_.rec[2]),"%s",NRTK_VERSION);
    snprintf(ropt_.ant[1],sizeof(ropt_.ant[1]),"ADVNULLANTENNA");
    for (int i=0;i<3;i++) ropt_.apppos[i]=pos[i];
    snprintf(ropt_.comment[0],sizeof(ropt_.comment[0]),"virtual reference station (network RTK)");
    snprintf(ropt_.comment[1],sizeof(ropt_.comment[1]),"null antenna at marker, ITRF %s",cfg.frame);
    struct {int rs; int idx;} m[]={{RNX_SYS_GPS,NRTK_GPS},{RNX_SYS_GAL,NRTK_GAL},
                                   {RNX_SYS_CMP,NRTK_BDS},{RNX_SYS_QZS,NRTK_QZS}};
    for (auto &k:m) {
        char list[128],*p,*q;
        int n=0;
        snprintf(list,sizeof(list),"%s",cfg.vrs.plan[k.idx].out);
        for (p=list;*p;p=q) {
            while (*p==' ') p++;
            for (q=p;*q&&*q!=' ';q++) ;
            if (*q) *q++='\0';
            if (!*p) continue;
            for (const char *t: {"C","L","D","S"}) {
                if (n<MAXOBSTYPE) snprintf(ropt_.tobs[k.rs][n++],4,"%s%s",t,p);
            }
        }
        ropt_.nobs[k.rs]=n;
    }
    rs_=new RtcmState;
    nrtk_rtcmout_init(&rs_->enc,staid,cfg.msm,pos,1);
}
VrsWriter::~VrsWriter() { close(); }

void VrsWriter::close()
{
    if (rnx_) fclose(rnx_);
    if (rtcm_) fclose(rtcm_);
    if (diag_) fclose(diag_);
    if (epo_) fclose(epo_);
    rnx_=rtcm_=diag_=epo_=nullptr;
    if (rs_) {nrtk_rtcmout_free(&rs_->enc); delete rs_; rs_=nullptr;}
}

void VrsWriter::write(const NetEpoch &ep, const VrsOutput &out)
{
    std::string ts=time_str(ep.t,1);
    char id[8];
    if (epo_) {
        fprintf(epo_,"%s,%d,%s,%d,%d,%d,%.4f,%.4f,%d,%.1f\n",ts.c_str(),out.mode,
                out.mode?cfg_.sta[out.mode==1?out.cell.master:out.cell.nearest].id:"",
                out.n,out.corr.nok,out.info.nlockreset,out.iono_activity,
                out.corr.trop_resid,out.degraded,ep.proc_ms);
    }
    if (diag_&&out.mode==1) {
        for (int i=0;i<out.corr.ns;i++) {
            const nrtk_vsat_t *s=out.corr.sat+i;
            satno2id(s->sat,id);
            fprintf(diag_,"%s,%s,%s,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.6f,%d,%d\n",
                    ts.c_str(),id,cfg_.sta[out.cell.master].id,s->isref,s->ok,s->nbl,
                    s->dd_iono,s->dd_nondisp,s->dd_trop_model,s->rms_iono,s->rms_nd,
                    s->coef_iono[0],s->coef_iono[1],s->coef_nd[0],s->coef_nd[1],
                    s->extrap,s->newfix);
        }
    }
    if (out.mode!=1||out.n<=0) return;   /* only VRS epochs go to the files */

    if (nep_==0&&rnx_) {
        ropt_.tstart=ep.t;
        outrnxobsh(rnx_,&ropt_,ep.nav.get());
    }
    if (rnx_) outrnxobsb(rnx_,&ropt_,out.obs,out.n,0);
    if (rtcm_) {
        uint8_t buf[32768];
        int n=0;
        if (nep_%10==0) n+=nrtk_rtcmout_station(&rs_->enc,buf+n,(int)sizeof(buf)-n);
        n+=nrtk_rtcmout_obs(&rs_->enc,ep.t,out.obs,out.n,buf+n,(int)sizeof(buf)-n);
        fwrite(buf,1,n,rtcm_);
    }
    nep_++;
}

} // namespace nrtk
