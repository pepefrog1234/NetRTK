/*------------------------------------------------------------------------------
 * nrtk/atmos.h : atmospheric error interpolation (spec ch.6)
 *
 * For every satellite the DD dispersive (ionosphere, L1) and non-dispersive
 * (troposphere + orbit) errors of the baselines master->aux_i are fitted with
 * a plane through the master in local north/east coordinates
 *      dC(M,A_i) = a*dN_i + b*dE_i,      dC(M,V) = a*dN_V + b*dE_V
 * n=2 auxiliary stations: linear interpolation model (LIM, exact),
 * n>=3: least-squares model (LSM) with residual RMS as quality indicator.
 * The non-dispersive part is reduced by the DD Saastamoinen+GMF model before
 * the fit and the model is restored with the height of the VRS.
 *----------------------------------------------------------------------------*/
#ifndef NRTK_ATMOS_H
#define NRTK_ATMOS_H
#include "nrtk/common.h"
#include "nrtk/prep.h"
#include "nrtk/netpos.h"
#ifdef __cplusplus
extern "C" {
#endif

enum { NRTK_INTERP_AUTO=0, NRTK_INTERP_LIM=1, NRTK_INTERP_LSM=2 };

typedef struct {
    int    method;              /* NRTK_INTERP_xxx */
    int    min_bl;              /* min baselines per satellite (2) */
    int    trop_model;          /* reduce non-dispersive part by the model */
    double iono_rate_thres;     /* ionosphere activity threshold (m/min) */
    double max_extrap_km;       /* max extrapolation distance (km) */
    int    min_sats;            /* fewer usable satellites -> degraded */
} nrtk_atmos_opt_t;

void nrtk_atmos_opt_default(nrtk_atmos_opt_t *opt);

/* one baseline of a cell as seen from the master */
typedef struct {
    const nrtk_blepoch_t *bl;   /* baseline epoch result */
    int    flip;                /* 1: bl is stored aux->master, values negated */
    int    aux;                 /* aux station index */
    double apos[3];             /* aux marker ECEF (m) */
} nrtk_cellbl_t;

typedef struct {
    uint8_t sat;
    int8_t  sysidx;
    uint8_t ok;                 /* correction available */
    uint8_t isref;              /* reference satellite */
    uint8_t nbl;                /* baselines used */
    uint8_t extrap;             /* VRS outside the hull of the used stations */
    uint8_t newfix;             /* a contributing ambiguity was (re)fixed */
    uint32_t blmask;            /* contributing baselines (bit i: cell baseline i) */
    double  dd_iono;            /* DD ionosphere master->VRS (L1, m) */
    double  dd_nondisp;         /* DD non-dispersive master->VRS (m) */
    double  dd_trop_model;      /* model part of dd_nondisp (m) */
    double  rms_iono,rms_nd;    /* LSM fit residual RMS (m), 0 for LIM */
    double  coef_iono[2];       /* plane coefficients north/east (m/km) */
    double  coef_nd[2];
} nrtk_vsat_t;

enum { NRTK_DEG_IONO=1, NRTK_DEG_FEWSAT=2, NRTK_DEG_EXTRAP=4, NRTK_DEG_NEAREST=8 };

typedef struct {
    gtime_t time;
    int     master;
    int     refsat[NRTK_NSYS];
    int     ns;
    nrtk_vsat_t sat[NRTK_BL_MAXSAT];
    int     nok;                /* satellites with corrections (incl. refs) */
    double  iono_activity;      /* max |d(DD iono)/dt| over the cell (m/min) */
    double  trop_resid;         /* RMS of reduced non-dispersive residuals (m) */
    int     degraded;           /* NRTK_DEG_xxx flags */
} nrtk_vcorr_t;

/* interpolate corrections to vpos. master: prepared master epoch (satellite
 * positions), bls: baselines master->aux of the cell. */
int nrtk_atmos_interp(const nrtk_atmos_opt_t *opt, const nrtk_staepoch_t *master,
                      const nrtk_cellbl_t *bls, int nbl, const double *vpos,
                      nrtk_vcorr_t *out);

/* ionosphere activity monitor: tracks DD ionosphere rates per baseline */
typedef struct nrtk_ionomon nrtk_ionomon_t;
nrtk_ionomon_t *nrtk_ionomon_create(double window_s);
void   nrtk_ionomon_free(nrtk_ionomon_t *m);
void   nrtk_ionomon_add(nrtk_ionomon_t *m, const nrtk_blepoch_t *bl);
double nrtk_ionomon_rate(const nrtk_ionomon_t *m, int sta_m, int sta_a); /* m/min */

/* plane fit helper (exposed for tests): y_i = a*n_i + b*e_i, n>=2.
 * returns 0 on failure, rms of residuals in *rms (0 for n=2). */
int nrtk_plane_fit(const double *n, const double *e, const double *y, int cnt,
                   double *a, double *b, double *rms);

#ifdef __cplusplus
}
#endif
#endif
