/*------------------------------------------------------------------------------
 * nrtk/prep.h : station epoch preparation
 *
 * Converts one epoch of raw observations of a reference station into the
 * "prepared" form used by all downstream modules:
 *   - signals reordered by the signal plan: slots 0..nproc-1 hold the
 *     processing frequencies of the system, further slots hold other signals
 *     that are carried to the VRS output
 *   - receiver antenna PCO/PCV and ARP eccentricity removed (marker, null
 *     antenna), phase in cycles, code in metres
 *   - receiver clock bias estimated from code with the known position
 *   - satellite positions at transmission time, geometric ranges to the
 *     tide-displaced marker, azimuth/elevation
 *----------------------------------------------------------------------------*/
#ifndef NRTK_PREP_H
#define NRTK_PREP_H
#include "nrtk/common.h"
#include "nrtk/antenna.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint8_t sat;
    int     sysidx;                 /* NRTK_xxx */
    int     nproc;                  /* processing slots present for this system */
    double  rs[6];                  /* satellite pos/vel at transmission (m,m/s) */
    double  dts[2];                 /* satellite clock bias/drift (s,s/s) */
    int     svh;
    double  azel[2];                /* azimuth/elevation at the station (rad) */
    double  r;                      /* geometric range incl. Sagnac (m) */
    double  rate;                   /* geometric range rate (m/s) */
    double  freq[NRTK_MAXSIG];      /* carrier frequency per slot (Hz) */
    double  dant[NRTK_MAXSIG];      /* removed antenna range correction (m) */
    obsd_t  obs;                    /* prepared observation (see above) */
} nrtk_satdata_t;

typedef struct {
    gtime_t time;                   /* receiver time tag (GPST) */
    int     sta;                    /* station index */
    double  rr[3];                  /* tide-displaced marker ECEF (m) */
    double  dtr;                    /* receiver clock bias (s) */
    int     ns;
    nrtk_satdata_t sat[NRTK_MAXSATEP];
} nrtk_staepoch_t;

typedef struct {
    double elmin;                   /* elevation mask (rad) */
    int    sysmask;                 /* enabled systems (1<<NRTK_xxx) */
    int    tide;                    /* apply solid earth tide */
    int    pcv;                     /* apply PCV in addition to PCO */
    nrtk_sigplan_t plan[NRTK_NSYS];
} nrtk_prep_opt_t;

void nrtk_prep_opt_default(nrtk_prep_opt_t *opt);

/* prepare one epoch. obs: raw observations of one station at one epoch.
 * ant may be NULL (no antenna model). Returns number of satellites. */
int nrtk_prep_epoch(const nrtk_prep_opt_t *opt, const nrtk_station_t *sta,
                    int ista, const nrtk_ant_t *ant, const obsd_t *obs, int n,
                    const nav_t *nav, nrtk_staepoch_t *ep);

/* find satellite in prepared epoch, NULL if absent */
const nrtk_satdata_t *nrtk_ep_sat(const nrtk_staepoch_t *ep, int sat);

/* phase (m) / code (m) of prepared slot, 0 if missing */
double nrtk_phase_m(const nrtk_satdata_t *s, int slot);
double nrtk_code_m (const nrtk_satdata_t *s, int slot);

#ifdef __cplusplus
}
#endif
#endif
