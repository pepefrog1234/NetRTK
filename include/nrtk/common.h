/*------------------------------------------------------------------------------
 * nrtk/common.h : libnrtk common definitions
 *
 * libnrtk is a network-RTK (virtual reference station) core built on top of
 * RTKLIB-EX. Modules exchange data only through RTKLIB public structures
 * (obsd_t, nav_t, rtk_t, rtcm_t) and the plain C structures defined here.
 *
 * Conventions used throughout the library
 *  - Single difference (SD) between stations X and Y of a quantity q:
 *        SD_XY(q) = q_Y - q_X          ("from X to Y")
 *  - Double difference between satellites s and reference r:
 *        DD_XY(q,s,r) = SD_XY(q,s) - SD_XY(q,r)
 *  - Baseline (M,A): M = master side (RTKLIB "base", rcv=2),
 *                    A = auxiliary side (RTKLIB "rover", rcv=1).
 *  - Ionospheric delays are expressed as the delay on the GPS L1 frequency
 *    (FREQL1, metres). The delay on frequency f is I*(FREQL1/f)^2.
 *  - Observations handed between modules are "prepared": antenna phase
 *    centre offsets/variations and the ARP eccentricity are removed, i.e.
 *    they refer to the station marker with a null antenna (see prep.h).
 *----------------------------------------------------------------------------*/
#ifndef NRTK_COMMON_H
#define NRTK_COMMON_H

#include <stdint.h>
#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NRTK_VERSION      "0.1.0"
#define NRTK_MAXSTA       160           /* max number of reference stations */
#define NRTK_MAXSIG       (NFREQ+NEXOBS) /* signal slots per satellite */
#define NRTK_NSYS         4             /* GPS, Galileo, BeiDou, QZSS */
#define NRTK_NPROC        3             /* max processing frequencies per system */
#define NRTK_MAXSATEP     MAXOBS        /* max satellites per station epoch */

/* system index used inside libnrtk */
enum { NRTK_GPS=0, NRTK_GAL=1, NRTK_BDS=2, NRTK_QZS=3 };

int         nrtk_sysidx(int sys);            /* SYS_xxx -> NRTK_xxx, -1: unsupported */
int         nrtk_idxsys(int idx);            /* NRTK_xxx -> SYS_xxx */
const char *nrtk_sysname(int idx);           /* "GPS","GAL","BDS","QZS" */
int         nrtk_sysmask(const char *list);  /* "GPS,GAL" -> mask of 1<<NRTK_xxx */

/* signal plan ---------------------------------------------------------------
 * For each system: which RINEX signals feed the processing frequencies
 * (slots 0..nproc-1, in priority order per frequency) and which signals are
 * carried through to the VRS output.                                         */
typedef struct {
    int    sys;                         /* SYS_xxx */
    int    nproc;                       /* processing frequencies (2 or 3) */
    char   proc[NRTK_NPROC][32];        /* e.g. "1C", "2W 2L 2X 2S", "5Q 5X 5I" */
    double pfreq[NRTK_NPROC];           /* carrier frequency of each proc slot (Hz) */
    char   out[96];                     /* output signals, e.g. "1C 2W 2L 5Q" */
} nrtk_sigplan_t;

void   nrtk_sigplan_default(nrtk_sigplan_t plan[NRTK_NSYS]);
double nrtk_code2freq(int sys, uint8_t code); /* carrier frequency (Hz), 0: n/a */
int    nrtk_sig_in_list(uint8_t code, const char *list);

/* reference station ---------------------------------------------------------*/
typedef struct {
    char   id[16];                      /* station id (4-9 chars) */
    double pos[3];                      /* marker ECEF at coord_epoch (m) */
    double vel[3];                      /* station velocity ECEF (m/yr) */
    double coord_epoch;                 /* coordinate epoch (decimal year) */
    char   frame[16];                   /* reference frame, e.g. ITRF2020 */
    char   antenna[64];                 /* IGS antenna type + radome */
    double antdel[3];                   /* ARP offset from marker E/N/U (m) */
    char   source[512];                 /* real-time stream (ntrip://,tcpcli://,file://,rinex://) */
    char   obsfile[512];                /* offline RINEX OBS file */
    int    rtcm_id;                     /* RTCM station id of the physical station */
    int    enabled;
} nrtk_station_t;

/* station marker position at time t (applies velocity), ECEF (m) */
void   nrtk_sta_pos(const nrtk_station_t *sta, gtime_t t, double *rr);
double nrtk_decimal_year(gtime_t t);

/* geometry helpers ----------------------------------------------------------*/
/* satellite position/clock at signal transmission for a receiver at rr that
 * receives at true time trx (light-time iteration, Sagnac included).
 * returns geometric range (m), <=0 on error. rs:{pos,vel} (m,m/s). */
double nrtk_satrange(gtime_t trx, const double *rr, int sat, const nav_t *nav,
                     double *rs, double *dts, double *rate, int *svh);
/* solid earth tide displaced position (ECEF, m) */
void   nrtk_tide_pos(gtime_t t, const double *r0, double *rr);
/* local north/east (m) of r relative to origin r0 */
void   nrtk_ne(const double *r0, const double *r, double *dn, double *de);
double nrtk_dist(const double *a, const double *b);

/* logging -------------------------------------------------------------------*/
enum { NRTK_LOG_ERROR=0, NRTK_LOG_WARN=1, NRTK_LOG_INFO=2, NRTK_LOG_DEBUG=3 };
void nrtk_log_level(int level);
void nrtk_log_json(int on);            /* structured (JSON lines) output */
void nrtk_log(int level, const char *event, const char *fmt, ...);
int  nrtk_log_enabled(int level);

#ifdef __cplusplus
}
#endif
#endif /* NRTK_COMMON_H */
