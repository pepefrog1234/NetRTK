/*------------------------------------------------------------------------------
 * nrtk/config.h : configuration (YAML) - loading and validation
 *
 * Example (spec ch.8):
 *   network:  {frame: ITRF2020, epoch: 2026.75, max_baseline_km: 70, antex: igs20.atx}
 *   stations: [{id: STA1, ecef: [x,y,z], antenna: "SEPCHOKE_B3E6  NONE",
 *               arp_height_m: 0.0, stream: ntrip://u:p@caster:2101/STA1,
 *               rinex: STA1.rnx}]
 *   solver:   {systems: [GPS,GAL,BDS,QZS], elevation_mask_deg: 10,
 *              ar_ratio: 3.0, iono_prior_ppm: 0.5}
 *   output:   {msm: 7, rate_hz: 1, vrs_rebuild_km: 5}
 * See conf/example.yaml for every key.
 *----------------------------------------------------------------------------*/
#ifndef NRTK_CONFIG_H
#define NRTK_CONFIG_H
#include "nrtk/common.h"
#include "nrtk/prep.h"
#include "nrtk/netpos.h"
#include "nrtk/atmos.h"
#include "nrtk/cellmgr.h"
#include "nrtk/vrsgen.h"
#ifdef __cplusplus
extern "C" {
#endif

#define NRTK_MAXJOB  32
#define NRTK_MAXUSER 64

typedef struct {
    char   name[64];            /* job name (output file prefix) */
    double pos[3];              /* VRS position ECEF (m), valid if has_pos */
    int    has_pos;
    char   leave_out[16];       /* station id: VRS at this station, station excluded */
    int    staid;               /* RTCM station id of the VRS (0: auto) */
} nrtk_vrsjob_t;

typedef struct {
    char user[32],pass[64];
    int  max_conn;
} nrtk_user_t;

/* mountpoints whose RTCM 1006 position is shifted to an output datum */
typedef struct {
    int  datum;                 /* index into nrtk_config_t.datum */
    char msm7[32],msm4[32],near[32];   /* empty: not offered */
} nrtk_dmount_t;

enum { NRTK_DATUM_INTERP_BARY=0, NRTK_DATUM_INTERP_NEAREST=1 };

typedef struct {
    /* network ---------------------------------------------------------------*/
    char   name[64];
    char   frame[32];
    double coord_epoch;
    char   antex[1024];
    int    nsta;
    nrtk_station_t sta[NRTK_MAXSTA];
    /* output datums (datum.h) ----------------------------------------------*/
    int    ndatum;
    nrtk_datum_def_t datum[NRTK_MAXDATUM];
    int    datum_interp;            /* NRTK_DATUM_INTERP_xxx */
    double datum_extrap_max_km;     /* beyond: nearest point's shift (20) */
    double datum_cover_km;          /* calibration point this close covers a station without coordinates (10) */
    char   datum_points_file[1024];
    int    ndpoint;
    nrtk_dpoint_t dpoint[NRTK_MAXDPOINT];
    /* processing -----------------------------------------------------------*/
    nrtk_prep_opt_t    prep;
    nrtk_netpos_opt_t  netpos;
    nrtk_atmos_opt_t   atmos;
    nrtk_cellmgr_opt_t cell;
    nrtk_vrs_opt_t     vrs;         /* template (position set per session) */
    /* output ---------------------------------------------------------------*/
    int    msm;                     /* 4 or 7 */
    double rate_hz;
    double vrs_rebuild_km;
    int    eph_out;                 /* pass ephemerides through (RTCM) */
    int    diag_csv;                /* write diagnostics CSV */
    /* offline --------------------------------------------------------------*/
    char   nav[4096];               /* navigation files */
    gtime_t ts,te;                  /* processing window ({0}: whole files) */
    double tint;                    /* epoch interval (s) */
    int    njob;
    nrtk_vrsjob_t job[NRTK_MAXJOB];
    char   outdir[1024];
    /* real-time server -----------------------------------------------------*/
    int    port;                    /* NTRIP caster port */
    char   mount_msm7[32],mount_msm4[32],mount_near[32];
    int    enable_near;
    int    ndmount;
    nrtk_dmount_t dmount[NRTK_MAXDATUM];
    double gga_timeout_s;           /* close if no valid GGA (15) */
    double nodata_close_s;          /* close if no VRS can be generated (10) */
    double sync_window_s;           /* epoch wait window (0.5) */
    double station_offline_s;       /* absent this long -> offline (10) */
    double health_fixrate;          /* station health: fix rate (0.9) */
    double health_latency_s;        /* station health: latency (2.0) */
    int    threads;                 /* baseline worker threads (0: auto) */
    int    allow_anonymous;
    int    max_conn_per_user;
    int    nuser;
    nrtk_user_t user[NRTK_MAXUSER];
    char   log_file[1024];
    int    log_json;
    double replay_speed;            /* rinex:// replay speed (1: real time) */
} nrtk_config_t;

void nrtk_config_default(nrtk_config_t *cfg);
/* load YAML file. returns 1 ok, 0 error (message in err) */
int  nrtk_config_load(const char *file, nrtk_config_t *cfg, char *err, int errlen);
/* semantic validation. returns number of errors (messages in err) */
int  nrtk_config_validate(const nrtk_config_t *cfg, char *err, int errlen);
int  nrtk_config_station(const nrtk_config_t *cfg, const char *id);
int  nrtk_config_datum(const nrtk_config_t *cfg, const char *name);   /* -1: unknown */
/* parse "YYYY-MM-DD hh:mm:ss" (GPST) */
int  nrtk_parse_time(const char *s, gtime_t *t);

#ifdef __cplusplus
}
#endif
#endif
