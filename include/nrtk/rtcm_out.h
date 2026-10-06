/*------------------------------------------------------------------------------
 * nrtk/rtcm_out.h : RTCM 3 output (MSM4/MSM7, 1006, 1033, ephemerides)
 *
 * MSM messages are split so that nsat*nsig <= 63 per message (RTKLIB-EX
 * drops the data of the 64th cell), with the multiple message bit set on all
 * but the last message of an epoch.
 * MSM and 1033 are encoded with RTKLIB-EX rtcm3e.c (one encoder state per
 * stream: lock times and phase-range integer offsets are stateful). 1006 is
 * encoded here because the RTKLIB encoder always flags a physical station and
 * ITRF year 0; a VRS must set DF141 (non-physical reference station) = 1.
 *----------------------------------------------------------------------------*/
#ifndef NRTK_RTCM_OUT_H
#define NRTK_RTCM_OUT_H
#include "nrtk/common.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    rtcm_t *rtcm;               /* RTKLIB encoder state */
    int     msm;                /* 4 or 7 */
    int     staid;              /* reference station id (0-4095) */
    double  pos[3];             /* ARP ECEF (m) */
    double  hgt;                /* antenna height above marker (m) */
    int     non_physical;       /* DF141 */
    int     itrf_year;          /* DF021 (0: not specified) */
    char    antdes[32];         /* antenna descriptor */
    char    rectype[32],recver[32],recsno[32];
    int     sysmask;            /* systems to encode (1<<NRTK_xxx) */
} nrtk_rtcmout_t;

int  nrtk_rtcmout_init(nrtk_rtcmout_t *o, int staid, int msm, const double *pos,
                       int non_physical);
void nrtk_rtcmout_free(nrtk_rtcmout_t *o);
/* change station id/position (VRS rebuild): resets lock times */
void nrtk_rtcmout_restart(nrtk_rtcmout_t *o, int staid, const double *pos);

/* 1006 + 1033. returns bytes written */
int nrtk_rtcmout_station(nrtk_rtcmout_t *o, uint8_t *buf, int size);
/* MSM messages for all systems in obs (multiple message bit set except last).
 * returns bytes written */
int nrtk_rtcmout_obs(nrtk_rtcmout_t *o, gtime_t t, const obsd_t *obs, int n,
                     uint8_t *buf, int size);
/* broadcast ephemeris message for one satellite (1019/1042/1044/1046) */
int nrtk_rtcmout_eph(nrtk_rtcmout_t *o, const nav_t *nav, int sat,
                     uint8_t *buf, int size);

/* encode 1006 into buf (exposed for tests). returns bytes written */
int nrtk_rtcm_encode_1006(int staid, const double *pos, double hgt,
                          int non_physical, int itrf_year, int sysmask,
                          uint8_t *buf, int size);

#ifdef __cplusplus
}
#endif
#endif
