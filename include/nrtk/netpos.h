/*------------------------------------------------------------------------------
 * nrtk/netpos.h : network solution - ambiguity fixing between reference
 *                 stations with known coordinates and residual extraction
 *
 * Two float/AR engines share the layers below (solver.engine):
 *   rtklib : RTKLIB-EX rtkpos() in fixed-position mode as the float filter
 *            (spec ch.5 primary design: SD ionosphere states, ZWD, slip
 *            handling), integer search/hold done here from outside
 *   native : the spec's back-up "own netpos filter": ionosphere-free phase
 *            with the fixed wide-lane removed, one float L1 ambiguity per
 *            satellite plus one relative ZWD per baseline (no code, no
 *            ionosphere states)
 * Shared layers
 *   1. extra-wide-lane (f2-f3) and wide-lane (f1-f2) fixing from smoothed
 *      Melbourne-Wubbena combinations,
 *   2. L1 ambiguity resolution with LAMBDA conditioned on the fixed wide-lanes,
 *      partial fixing (high-elevation satellites first), ratio test,
 *   3. fix-and-hold of the fixed ambiguities into the RTKLIB filter state,
 *   4. MW / geometry-free slip detection injected through LLI,
 *   5. extraction of ambiguity-fixed DD phase residuals and their split into
 *      the dispersive (ionosphere) and non-dispersive part.
 * The RTKLIB sources are used unmodified.
 *----------------------------------------------------------------------------*/
#ifndef NRTK_NETPOS_H
#define NRTK_NETPOS_H
#include "nrtk/common.h"
#include "nrtk/prep.h"
#ifdef __cplusplus
extern "C" {
#endif

enum { NRTK_ENGINE_NATIVE=0, NRTK_ENGINE_RTKLIB=1 };

typedef struct {
    int    engine;             /* NRTK_ENGINE_xxx */
    int    sysmask;            /* systems used (1<<NRTK_xxx) */
    int    nfreq;              /* processing frequencies: 2 or 3 (EWL) */
    double elmin;              /* elevation mask (rad) */
    double elmin_ar;           /* elevation mask for fixing new ambiguities (rad) */
    double ar_ratio;           /* LAMBDA ratio threshold */
    int    min_fix;            /* min DD ambiguities for an initial fix */
    int    min_lock;           /* min continuous epochs before AR */
    int    max_drop;           /* max satellites dropped in partial AR */
    double iono_prior_ppm;     /* SD ionosphere prior std (ppm of baseline) */
    double iono_prior_min;     /* lower bound of the prior std (m) */
    double iono_prn;           /* iono process noise (m/sqrt(s)) per 10 km */
    double trop_prn;           /* ZWD process noise (m/sqrt(s)) */
    double trop_std;           /* initial ZWD std (m) */
    int    trop_grad;          /* estimate gradients (est-ztdgrad) */
    double err_phase;          /* phase error (m), elevation model a+b/sin(el) */
    double code_ratio;         /* code/phase error ratio */
    double wl_min_epochs;      /* min MW samples before WL fix */
    double wl_min_time;        /* min MW averaging time before WL fix (s) */
    double wl_max_frac;        /* max |fractional part| to fix WL (cycles) */
    double wl_max_std;         /* max std of mean MW (cycles) */
    double mw_slip;            /* MW slip threshold (cycles) */
    double gf_slip;            /* geometry-free slip threshold (m) */
    double hold_var;           /* fix-and-hold pseudo-measurement variance (cyc^2) */
    double max_res;            /* max |DD L1/L2 residual - prediction| for a fixed sat (m) */
    int    tide;               /* solid earth tide */
    int    closure;            /* network: triangle closure check of the fixes */
    int    trace;              /* RTKLIB trace level (0: off) */
} nrtk_netpos_opt_t;

void nrtk_netpos_opt_default(nrtk_netpos_opt_t *opt);

/* per satellite output of one baseline epoch */
typedef struct {
    uint8_t sat;
    int8_t  sysidx;
    uint8_t fixed;              /* f1/f2 ambiguities fixed (held) */
    uint8_t wlfixed;            /* wide-lane fixed */
    uint8_t ewlfixed;           /* extra-wide-lane fixed */
    uint8_t isref;              /* reference satellite of its system */
    uint8_t newfix;             /* fixed (or re-fixed) in this epoch */
    uint8_t slip;               /* slip detected this epoch */
    float   azel_m[2],azel_a[2];/* az/el at master and aux (rad) */
    double  dd_res[NRTK_NPROC]; /* fixed DD phase residual per proc freq (m) */
    double  dd_code[NRTK_NPROC];/* DD code residual (m) */
    double  dd_iono;            /* DD ionospheric delay on L1 (m) */
    double  dd_nondisp;         /* DD non-dispersive delay (m) */
    double  dd_trop_model;      /* DD model troposphere (m), for diagnostics */
    int     nfixep;             /* consecutive fixed epochs */
} nrtk_blsat_t;

#define NRTK_BL_MAXSAT 96

typedef struct {
    gtime_t time;
    int     sta_m,sta_a;        /* station indices: master side / aux side */
    int     valid;              /* epoch processed */
    int     stat;               /* 0:none 1:float 2:fixed */
    double  ratio;              /* last LAMBDA ratio */
    int     nsat;               /* common satellites above mask */
    int     nfixed;             /* satellites with fixed ambiguities */
    int     nwl;                /* satellites with fixed WL */
    int     nslip;              /* slips in this epoch */
    double  ztd[2];             /* estimated ZTD master/aux (m) */
    double  dt_age;             /* time difference between stations (s) */
    double  length;             /* baseline length (m) */
    int     refsat[NRTK_NSYS];  /* reference satellite per system (0: none) */
    int     ns;
    nrtk_blsat_t sat[NRTK_BL_MAXSAT];
} nrtk_blepoch_t;

typedef struct nrtk_baseline nrtk_baseline_t;

nrtk_baseline_t *nrtk_bl_create(const nrtk_netpos_opt_t *opt,
                                const nrtk_station_t *m, int ista_m,
                                const nrtk_station_t *a, int ista_a);
void nrtk_bl_free(nrtk_baseline_t *bl);
void nrtk_bl_reset(nrtk_baseline_t *bl);

/* process one synchronized epoch of the two stations (prepared epochs).
 * returns 1 if the baseline produced output. */
int nrtk_bl_process(nrtk_baseline_t *bl, const nrtk_staepoch_t *em,
                    const nrtk_staepoch_t *ea, const nav_t *nav,
                    nrtk_blepoch_t *out);

/* external integrity information (network closure): drop the integer fix of
 * satellite sat (and its wide-lane) so it is resolved again */
void nrtk_bl_unfix(nrtk_baseline_t *bl, int sat);

/* cumulative statistics */
typedef struct {
    long   nepoch;              /* processed epochs */
    long   nfixep;              /* fixed epochs */
    double ttff;                /* time to first fix (s), <0: not fixed */
    long   nrefix;              /* fix losses followed by refix */
    long   nslip;
    double fixrate_60s;         /* fix rate over the last 60 s */
} nrtk_blstat_t;

void nrtk_bl_stat(const nrtk_baseline_t *bl, nrtk_blstat_t *st);

/* helpers shared with atmos/tests --------------------------------------------*/
/* split DD phase residuals on two frequencies into L1 ionosphere and
 * non-dispersive delay: d_a = D - I*(f_L1/f_a)^2, d_b = D - I*(f_L1/f_b)^2 */
void nrtk_split_iono(double fa, double fb, double da, double db,
                     double *iono_l1, double *nondisp);

#ifdef __cplusplus
}
#endif
#endif
