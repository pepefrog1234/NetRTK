/*------------------------------------------------------------------------------
 * config.cpp : YAML configuration loading and validation (yaml-cpp)
 *----------------------------------------------------------------------------*/
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>
#include "nrtk/config.h"

namespace {

void copy(char *dst, size_t n, const std::string &s)
{
    snprintf(dst,n,"%s",s.c_str());
}
template<class T> bool get(const YAML::Node &n, const char *key, T &v)
{
    if (!n||!n.IsMap()||!n[key]) return false;
    v=n[key].as<T>();
    return true;
}
bool get_str(const YAML::Node &n, const char *key, char *dst, size_t len)
{
    std::string s;
    if (!get(n,key,s)) return false;
    copy(dst,len,s);
    return true;
}
bool get_bool(const YAML::Node &n, const char *key, int &v)
{
    bool b;
    if (!get(n,key,b)) return false;
    v=b?1:0;
    return true;
}
bool get_deg(const YAML::Node &n, const char *key, double &rad)
{
    double d;
    if (!get(n,key,d)) return false;
    rad=d*D2R;
    return true;
}
std::string join_list(const YAML::Node &n)
{
    std::string s;
    if (n.IsSequence()) {
        for (auto v:n) s+=v.as<std::string>()+" ";
    }
    else s=n.as<std::string>();
    return s;
}
bool get_pos(const YAML::Node &n, double *ecef)
{
    if (n["ecef"]) {
        for (int i=0;i<3;i++) ecef[i]=n["ecef"][i].as<double>();
        return true;
    }
    if (n["llh"]) { /* deg, deg, m (ellipsoidal) */
        double pos[3];
        pos[0]=n["llh"][0].as<double>()*D2R;
        pos[1]=n["llh"][1].as<double>()*D2R;
        pos[2]=n["llh"][2].as<double>();
        pos2ecef(pos,ecef);
        return true;
    }
    return false;
}
/* official coordinates of a marker in the output datums:
 *   datum: {TWD97_2010: {ecef: [...], source: control-tie, observed: 2026-09-20, sigma_m: 0.015}} */
bool get_datums(const YAML::Node &n, const nrtk_config_t *c, const char *who,
                nrtk_stadatum_t *out, char *err, int errlen)
{
    YAML::Node dn=n["datum"];
    if (!dn) return true;
    if (!dn.IsMap()) {snprintf(err,errlen,"%s: datum must be a map of datum names",who); return false;}
    for (auto it:dn) {
        std::string name=it.first.as<std::string>();
        int k=nrtk_config_datum(c,name.c_str());
        if (k<0) {
            snprintf(err,errlen,"%s: datum %s is not in network.output_datums",who,name.c_str());
            return false;
        }
        nrtk_stadatum_t *d=out+k;
        memset(d,0,sizeof(*d));
        if (!get_pos(it.second,d->pos)) {
            snprintf(err,errlen,"%s: datum %s needs ecef or llh",who,name.c_str());
            return false;
        }
        get(it.second,"sigma_m",d->sigma);
        get_str(it.second,"source",d->source,sizeof(d->source));
        get_str(it.second,"observed",d->observed,sizeof(d->observed));
        d->valid=1;
    }
    return true;
}
/* calibration points file: a sequence (or {points: [...]}) of
 *   {id, ecef|llh, epoch, velocity, datum: {...}} */
bool load_dpoints(const char *file, nrtk_config_t *c, char *err, int errlen)
{
    YAML::Node y=YAML::LoadFile(file);
    YAML::Node seq=y.IsSequence()?y:y["points"];
    if (!seq||!seq.IsSequence()) {
        snprintf(err,errlen,"%s: expected a sequence of points",file); return false;
    }
    c->ndpoint=0;
    for (auto n:seq) {
        if (c->ndpoint>=NRTK_MAXDPOINT) {snprintf(err,errlen,"%s: too many points",file); return false;}
        nrtk_dpoint_t *p=c->dpoint+c->ndpoint;
        memset(p,0,sizeof(*p));
        if (!get_str(n,"id",p->id,sizeof(p->id))) {snprintf(err,errlen,"%s: point without id",file); return false;}
        if (!get_pos(n,p->pos)) {
            snprintf(err,errlen,"%s: point %s: ecef or llh (ITRF) required",file,p->id); return false;
        }
        if (n["velocity"]) for (int i=0;i<3;i++) p->vel[i]=n["velocity"][i].as<double>();
        p->coord_epoch=c->coord_epoch;
        get(n,"epoch",p->coord_epoch);
        std::string who=std::string(file)+": point "+p->id;
        if (!get_datums(n,c,who.c_str(),p->datum,err,errlen)) return false;
        c->ndpoint++;
    }
    return true;
}

} // namespace

extern "C" int nrtk_parse_time(const char *s, gtime_t *t)
{
    double ep[6]={0};
    if (sscanf(s,"%lf-%lf-%lf %lf:%lf:%lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<3&&
        sscanf(s,"%lf/%lf/%lf %lf:%lf:%lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<3) {
        return 0;
    }
    *t=epoch2time(ep);
    return 1;
}

extern "C" void nrtk_config_default(nrtk_config_t *c)
{
    memset(c,0,sizeof(*c));
    snprintf(c->name,sizeof(c->name),"NRTK");
    snprintf(c->frame,sizeof(c->frame),"ITRF2020");
    nrtk_prep_opt_default(&c->prep);
    nrtk_netpos_opt_default(&c->netpos);
    nrtk_atmos_opt_default(&c->atmos);
    nrtk_cellmgr_opt_default(&c->cell);
    nrtk_vrs_opt_default(&c->vrs);
    c->msm=7;
    c->rate_hz=1.0;
    c->vrs_rebuild_km=5.0;
    c->diag_csv=1;
    c->tint=1.0;
    snprintf(c->outdir,sizeof(c->outdir),"out");
    c->port=2101;
    snprintf(c->mount_msm7,sizeof(c->mount_msm7),"VRS_MSM7");
    snprintf(c->mount_msm4,sizeof(c->mount_msm4),"VRS_MSM4");
    snprintf(c->mount_near,sizeof(c->mount_near),"NEAR");
    c->enable_near=1;
    c->gga_timeout_s=15.0;
    c->nodata_close_s=10.0;
    c->sync_window_s=0.5;
    c->station_offline_s=10.0;
    c->health_fixrate=0.9;
    c->health_latency_s=2.0;
    c->allow_anonymous=1;
    c->max_conn_per_user=4;
    c->replay_speed=1.0;
    c->datum_interp=NRTK_DATUM_INTERP_BARY;
    c->datum_extrap_max_km=20.0;
    c->datum_cover_km=10.0;
}

extern "C" int nrtk_config_load(const char *file, nrtk_config_t *c, char *err,
                                int errlen)
{
    try {
        YAML::Node y=YAML::LoadFile(file);
        std::string s;

        /* network */
        YAML::Node nw=y["network"];
        get_str(nw,"name",c->name,sizeof(c->name));
        get_str(nw,"frame",c->frame,sizeof(c->frame));
        get(nw,"epoch",c->coord_epoch);
        if (get(nw,"max_baseline_km",c->cell.max_baseline_km)) {
            c->cell.max_dist_km=c->cell.max_baseline_km;
        }
        get(nw,"max_distance_km",c->cell.max_dist_km);
        get(nw,"lsm_radius_km",c->cell.lsm_radius_km);
        get(nw,"max_aux",c->cell.max_aux);
        get_str(nw,"antex",c->antex,sizeof(c->antex));

        /* output datums: [{name, label, version, frame, epoch}] */
        c->ndatum=0;
        if (nw&&nw["output_datums"]) {
            for (auto n:nw["output_datums"]) {
                if (c->ndatum>=NRTK_MAXDATUM) {snprintf(err,errlen,"too many output_datums (max %d)",NRTK_MAXDATUM); return 0;}
                nrtk_datum_def_t *d=c->datum+c->ndatum;
                memset(d,0,sizeof(*d));
                if (n.IsScalar()) copy(d->name,sizeof(d->name),n.as<std::string>());
                else {
                    get_str(n,"name",d->name,sizeof(d->name));
                    get_str(n,"label",d->label,sizeof(d->label));
                    get_str(n,"version",d->version,sizeof(d->version));
                    get_str(n,"frame",d->frame,sizeof(d->frame));
                    get(n,"epoch",d->epoch);
                }
                if (!d->name[0]) {snprintf(err,errlen,"output_datums: name required"); return 0;}
                for (const char *p=d->name;*p;p++) {
                    if (!isalnum((unsigned char)*p)&&*p!='_') {
                        snprintf(err,errlen,"output_datums: name %s: letters, digits and _ only",d->name);
                        return 0;
                    }
                }
                if (nrtk_config_datum(c,d->name)>=0) {snprintf(err,errlen,"output_datums: duplicate %s",d->name); return 0;}
                /* defaults of the theory chain for the Taiwan datums */
                if (!d->frame[0]) {
                    if (!strcmp(d->name,"TWD97_2010")) {strcpy(d->frame,"ITRF94"); d->epoch=2010.0;}
                    else if (!strcmp(d->name,"TWD97_2020")) {strcpy(d->frame,"ITRF2014"); d->epoch=2020.0;}
                    else if (!strcmp(d->name,"TWD97")) {strcpy(d->frame,"ITRF94"); d->epoch=1997.0;}
                }
                if (!d->label[0]) copy(d->label,sizeof(d->label),d->name);
                c->ndatum++;
            }
        }
        if (nw&&nw["datum_interp"]) {
            s=nw["datum_interp"].as<std::string>();
            if (s=="barycentric") c->datum_interp=NRTK_DATUM_INTERP_BARY;
            else if (s=="nearest") c->datum_interp=NRTK_DATUM_INTERP_NEAREST;
            else {snprintf(err,errlen,"network.datum_interp must be barycentric or nearest"); return 0;}
        }
        get(nw,"datum_extrap_max_km",c->datum_extrap_max_km);
        get(nw,"datum_cover_km",c->datum_cover_km);
        get_str(nw,"datum_points_file",c->datum_points_file,sizeof(c->datum_points_file));

        /* stations */
        c->nsta=0;
        if (y["stations"]) {
            for (auto n:y["stations"]) {
                if (c->nsta>=NRTK_MAXSTA) {snprintf(err,errlen,"too many stations"); return 0;}
                nrtk_station_t *st=c->sta+c->nsta;
                memset(st,0,sizeof(*st));
                st->enabled=1;
                if (!get_str(n,"id",st->id,sizeof(st->id))) {
                    snprintf(err,errlen,"station without id"); return 0;
                }
                if (!get_pos(n,st->pos)) {
                    snprintf(err,errlen,"station %s: ecef or llh required",st->id); return 0;
                }
                if (n["velocity"]) for (int i=0;i<3;i++) st->vel[i]=n["velocity"][i].as<double>();
                st->coord_epoch=c->coord_epoch;
                get(n,"epoch",st->coord_epoch);
                copy(st->frame,sizeof(st->frame),c->frame);
                get_str(n,"antenna",st->antenna,sizeof(st->antenna));
                get(n,"arp_height_m",st->antdel[2]);
                if (n["arp_enu_m"]) for (int i=0;i<3;i++) st->antdel[i]=n["arp_enu_m"][i].as<double>();
                get_str(n,"stream",st->source,sizeof(st->source));
                get_str(n,"rinex",st->obsfile,sizeof(st->obsfile));
                get(n,"rtcm_id",st->rtcm_id);
                get_bool(n,"enabled",st->enabled);
                std::string who=std::string("station ")+st->id;
                if (!get_datums(n,c,who.c_str(),st->datum,err,errlen)) return 0;
                c->nsta++;
            }
        }
        if (c->datum_points_file[0]&&!load_dpoints(c->datum_points_file,c,err,errlen)) return 0;
        /* solver */
        YAML::Node sv=y["solver"];
        if (sv&&sv["systems"]) {
            s=join_list(sv["systems"]);
            c->netpos.sysmask=c->prep.sysmask=c->vrs.sysmask=nrtk_sysmask(s.c_str());
        }
        if (get_deg(sv,"elevation_mask_deg",c->netpos.elmin)) c->vrs.elmin=c->netpos.elmin;
        get_deg(sv,"ar_elevation_mask_deg",c->netpos.elmin_ar);
        if (sv&&sv["engine"]) {
            s=sv["engine"].as<std::string>();
            if (s=="rtklib") c->netpos.engine=NRTK_ENGINE_RTKLIB;
            else if (s=="native") c->netpos.engine=NRTK_ENGINE_NATIVE;
            else {snprintf(err,errlen,"solver.engine must be native or rtklib"); return 0;}
        }
        get(sv,"frequencies",c->netpos.nfreq);
        get(sv,"ar_ratio",c->netpos.ar_ratio);
        get(sv,"min_fix_sats",c->netpos.min_fix);
        get(sv,"min_lock",c->netpos.min_lock);
        get(sv,"max_drop",c->netpos.max_drop);
        get(sv,"iono_prior_ppm",c->netpos.iono_prior_ppm);
        get(sv,"iono_prior_min_m",c->netpos.iono_prior_min);
        get(sv,"iono_process_noise",c->netpos.iono_prn);
        get(sv,"trop_process_noise",c->netpos.trop_prn);
        get(sv,"trop_initial_std_m",c->netpos.trop_std);
        if (sv&&sv["troposphere"]) {
            s=sv["troposphere"].as<std::string>();
            c->netpos.trop_grad=(s=="est-ztdgrad");
        }
        get(sv,"phase_error_m",c->netpos.err_phase);
        get(sv,"code_phase_ratio",c->netpos.code_ratio);
        get(sv,"wl_min_epochs",c->netpos.wl_min_epochs);
        get(sv,"wl_min_time_s",c->netpos.wl_min_time);
        get(sv,"wl_max_frac",c->netpos.wl_max_frac);
        get(sv,"mw_slip_cycles",c->netpos.mw_slip);
        get(sv,"gf_slip_m",c->netpos.gf_slip);
        get(sv,"hold_variance",c->netpos.hold_var);
        get(sv,"max_fix_residual_m",c->netpos.max_res);
        get_bool(sv,"tide",c->netpos.tide);
        get_bool(sv,"closure_check",c->netpos.closure);
        c->prep.tide=c->vrs.tide=c->netpos.tide;
        get(sv,"trace",c->netpos.trace);
        get_bool(sv,"pcv",c->prep.pcv);

        /* signal plan overrides: signals: {GPS: {proc: ["1C","2W 2L","5Q"], out: "1C 2W"}} */
        if (sv&&sv["signals"]) {
            for (auto it:sv["signals"]) {
                int idx=-1, mask=nrtk_sysmask(it.first.as<std::string>().c_str());
                for (int i=0;i<NRTK_NSYS;i++) if (mask==(1<<i)) idx=i;
                if (idx<0) continue;
                nrtk_sigplan_t *p=c->prep.plan+idx;
                if (it.second["proc"]) {
                    p->nproc=0;
                    for (auto f:it.second["proc"]) {
                        if (p->nproc>=NRTK_NPROC) break;
                        copy(p->proc[p->nproc++],sizeof(p->proc[0]),f.as<std::string>());
                    }
                }
                if (it.second["out"]) copy(p->out,sizeof(p->out),join_list(it.second["out"]));
            }
            memcpy(c->vrs.plan,c->prep.plan,sizeof(c->vrs.plan));
        }
        /* interpolation */
        YAML::Node ip=y["interp"];
        if (ip&&ip["method"]) {
            s=ip["method"].as<std::string>();
            c->atmos.method=s=="lim"?NRTK_INTERP_LIM:(s=="lsm"?NRTK_INTERP_LSM:NRTK_INTERP_AUTO);
        }
        get(ip,"iono_rate_threshold_m_per_min",c->atmos.iono_rate_thres);
        get(ip,"min_baselines",c->atmos.min_bl);
        get(ip,"min_sats",c->atmos.min_sats);
        get_bool(ip,"trop_model",c->atmos.trop_model);
        get_bool(ip,"height_scale",c->atmos.height_scale);
        get(ip,"max_extrapolation_km",c->atmos.max_extrap_km);

        /* output */
        YAML::Node out=y["output"];
        get(out,"msm",c->msm);
        get(out,"rate_hz",c->rate_hz);
        get(out,"vrs_rebuild_km",c->vrs_rebuild_km);
        get_bool(out,"ephemeris",c->eph_out);
        get_bool(out,"diagnostics",c->diag_csv);
        get_bool(out,"pass_unfixed",c->vrs.pass_unfixed);
        get_deg(out,"elevation_mask_deg",c->vrs.elmin);
        get(out,"lock_reset_jump_m",c->vrs.jump_thres);

        /* offline */
        YAML::Node of=y["offline"];
        if (of&&of["nav"]) copy(c->nav,sizeof(c->nav),join_list(of["nav"]));
        if (get(of,"start",s)&&!nrtk_parse_time(s.c_str(),&c->ts)) {
            snprintf(err,errlen,"bad offline.start: %s",s.c_str()); return 0;
        }
        if (get(of,"end",s)&&!nrtk_parse_time(s.c_str(),&c->te)) {
            snprintf(err,errlen,"bad offline.end: %s",s.c_str()); return 0;
        }
        double dur=0.0;
        if (get(of,"duration_s",dur)&&c->ts.time) c->te=timeadd(c->ts,dur);
        get(of,"interval_s",c->tint);
        get_str(of,"out_dir",c->outdir,sizeof(c->outdir));
        if (of&&of["vrs"]) {
            c->njob=0;
            for (auto j:of["vrs"]) {
                if (c->njob>=NRTK_MAXJOB) break;
                nrtk_vrsjob_t *jb=c->job+c->njob;
                memset(jb,0,sizeof(*jb));
                get_str(j,"name",jb->name,sizeof(jb->name));
                get_str(j,"leave_out",jb->leave_out,sizeof(jb->leave_out));
                jb->has_pos=get_pos(j,jb->pos);
                get(j,"rtcm_id",jb->staid);
                if (!jb->name[0]) snprintf(jb->name,sizeof(jb->name),"%s%s",
                    jb->leave_out[0]?"loo_":"vrs",jb->leave_out[0]?jb->leave_out:"");
                c->njob++;
            }
        }
        /* server */
        YAML::Node se=y["server"];
        get(se,"port",c->port);
        get_str(se,"mount_msm7",c->mount_msm7,sizeof(c->mount_msm7));
        get_str(se,"mount_msm4",c->mount_msm4,sizeof(c->mount_msm4));
        get_str(se,"mount_nearest",c->mount_near,sizeof(c->mount_near));
        get_bool(se,"nearest_mountpoint",c->enable_near);
        /* datum_mounts: [{datum: TWD97_2010, msm7: VRS_MSM7_TWD97_2010, msm4: ..., nearest: ...}] */
        c->ndmount=0;
        if (se&&se["datum_mounts"]) {
            for (auto n:se["datum_mounts"]) {
                if (c->ndmount>=NRTK_MAXDATUM) break;
                nrtk_dmount_t *m=c->dmount+c->ndmount;
                memset(m,0,sizeof(*m));
                if (!get(n,"datum",s)) {snprintf(err,errlen,"server.datum_mounts: datum required"); return 0;}
                if ((m->datum=nrtk_config_datum(c,s.c_str()))<0) {
                    snprintf(err,errlen,"server.datum_mounts: datum %s is not in network.output_datums",s.c_str());
                    return 0;
                }
                /* defaults: <mount>_<datum name> */
                snprintf(m->msm7,sizeof(m->msm7),"%s_%s",c->mount_msm7,s.c_str());
                snprintf(m->msm4,sizeof(m->msm4),"%s_%s",c->mount_msm4,s.c_str());
                snprintf(m->near,sizeof(m->near),"%s_%s",c->mount_near,s.c_str());
                get_str(n,"msm7",m->msm7,sizeof(m->msm7));
                get_str(n,"msm4",m->msm4,sizeof(m->msm4));
                get_str(n,"nearest",m->near,sizeof(m->near));
                c->ndmount++;
            }
        }
        get(se,"gga_timeout_s",c->gga_timeout_s);
        get(se,"no_data_close_s",c->nodata_close_s);
        get(se,"sync_window_s",c->sync_window_s);
        get(se,"station_offline_s",c->station_offline_s);
        get(se,"health_fixrate",c->health_fixrate);
        get(se,"health_latency_s",c->health_latency_s);
        get(se,"threads",c->threads);
        get_bool(se,"allow_anonymous",c->allow_anonymous);
        get(se,"max_connections_per_user",c->max_conn_per_user);
        get_str(se,"log_file",c->log_file,sizeof(c->log_file));
        get_bool(se,"log_json",c->log_json);
        get(se,"replay_speed",c->replay_speed);
        if (se&&se["users"]) {
            c->nuser=0;
            for (auto u:se["users"]) {
                if (c->nuser>=NRTK_MAXUSER) break;
                nrtk_user_t *us=c->user+c->nuser++;
                memset(us,0,sizeof(*us));
                get_str(u,"user",us->user,sizeof(us->user));
                get_str(u,"password",us->pass,sizeof(us->pass));
                us->max_conn=c->max_conn_per_user;
                get(u,"max_connections",us->max_conn);
            }
        }
    }
    catch (const std::exception &e) {
        snprintf(err,errlen,"%s: %s",file,e.what());
        return 0;
    }
    return 1;
}

extern "C" int nrtk_config_station(const nrtk_config_t *c, const char *id)
{
    for (int i=0;i<c->nsta;i++) if (!strcmp(c->sta[i].id,id)) return i;
    return -1;
}
extern "C" int nrtk_config_datum(const nrtk_config_t *c, const char *name)
{
    for (int i=0;i<c->ndatum;i++) if (!strcmp(c->datum[i].name,name)) return i;
    return -1;
}

extern "C" int nrtk_config_validate(const nrtk_config_t *c, char *err, int errlen)
{
    std::string msg;
    int nerr=0;
    auto add=[&](const std::string &m){ if (nerr++) msg+="; "; msg+=m; };

    if (c->nsta<1) add("no stations");
    for (int i=0;i<c->nsta;i++) {
        const nrtk_station_t *s=c->sta+i;
        double r=norm(s->pos,3);
        if (r<6.3E6||r>6.4E6) add(std::string("station ")+s->id+": position not on earth surface");
        for (int j=0;j<i;j++) {
            if (!strcmp(c->sta[j].id,s->id)) add(std::string("duplicate station id ")+s->id);
            else if (nrtk_dist(c->sta[j].pos,s->pos)<1.0) {
                add(std::string("stations ")+c->sta[j].id+" and "+s->id+" coincide");
            }
        }
        if (s->rtcm_id<0||s->rtcm_id>4095) add(std::string("station ")+s->id+": rtcm_id out of range");
        /* official datum coordinates: within 5 m of ITRF (a wrong version, unit or
           marker would be metres or kilometres off) */
        for (int k=0;k<c->ndatum;k++) {
            if (!s->datum[k].valid) continue;
            double dd=nrtk_dist(s->datum[k].pos,s->pos);
            if (dd>5.0) {
                char b[160];
                snprintf(b,sizeof(b),"station %s: %s coordinates %.1f m from the ITRF position",
                         s->id,c->datum[k].name,dd);
                add(b);
            }
        }
    }
    for (int i=0;i<c->ndpoint;i++) {
        const nrtk_dpoint_t *p=c->dpoint+i;
        double r=norm(p->pos,3);
        if (r<6.3E6||r>6.4E6) add(std::string("datum point ")+p->id+": position not on earth surface");
        if (nrtk_config_station(c,p->id)>=0) add(std::string("datum point ")+p->id+": same id as a station");
        for (int j=0;j<i;j++) if (!strcmp(c->dpoint[j].id,p->id)) add(std::string("duplicate datum point ")+p->id);
        int nd=0;
        for (int k=0;k<c->ndatum;k++) {
            if (!p->datum[k].valid) continue;
            nd++;
            if (nrtk_dist(p->datum[k].pos,p->pos)>5.0) {
                add(std::string("datum point ")+p->id+": "+c->datum[k].name+" coordinates far from ITRF");
            }
        }
        if (!nd) add(std::string("datum point ")+p->id+": no datum coordinates");
    }
    if (c->datum_extrap_max_km<0.0) add("network.datum_extrap_max_km must be >= 0");
    if (c->datum_cover_km<0.0) add("network.datum_cover_km must be >= 0");
    for (int i=0;i<c->ndmount;i++) for (int j=0;j<i;j++) {
        if (c->dmount[i].datum==c->dmount[j].datum) {
            add(std::string("server.datum_mounts: datum ")+c->datum[c->dmount[i].datum].name+" listed twice");
        }
    }
    if (c->msm!=4&&c->msm!=7) add("output.msm must be 4 or 7");
    if (c->netpos.ar_ratio<1.0) add("solver.ar_ratio must be >= 1");
    if (c->netpos.elmin<0.0||c->netpos.elmin>40.0*D2R) add("solver.elevation_mask_deg out of range");
    if (!c->netpos.sysmask) add("solver.systems empty");
    if (c->netpos.nfreq<2||c->netpos.nfreq>3) add("solver.frequencies must be 2 or 3");
    if (c->cell.max_baseline_km<=0.0) add("network.max_baseline_km must be > 0");
    for (int i=0;i<c->njob;i++) {
        const nrtk_vrsjob_t *j=c->job+i;
        if (j->leave_out[0]&&nrtk_config_station(c,j->leave_out)<0) {
            add(std::string("offline.vrs ")+j->name+": unknown station "+j->leave_out);
        }
        if (!j->leave_out[0]&&!j->has_pos) add(std::string("offline.vrs ")+j->name+": position required");
    }
    if (err&&errlen>0) snprintf(err,errlen,"%s",msg.c_str());
    return nerr;
}
