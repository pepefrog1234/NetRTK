/*------------------------------------------------------------------------------
 * nrtk/vrsgen.h : virtual observation generator (spec ch.7)
 *
 *   L_V,j = L_M,j + (rho_V - rho_M) + dD_MV - (f_L1/f_j)^2 dI_MV   (metres)
 *   P_V,j = P_M,j + (rho_V - rho_M) + dD_MV + (f_L1/f_j)^2 dI_MV
 *
 * L_M/P_M are the master's prepared (null antenna, marker) observations, so
 * the VRS refers to a null antenna at the VRS marker (ADVNULLANTENNA). The
 * master's receiver clock, ambiguities and slips are kept. The single
 * difference corrections d(.)_MV are formed from the interpolated DD values
 * with the reference satellite set to zero; when the reference changes the
 * per-system offset is carried over so that the corrections stay continuous.
 *----------------------------------------------------------------------------*/
#ifndef NRTK_VRSGEN_H
#define NRTK_VRSGEN_H
#include "nrtk/common.h"
#include "nrtk/prep.h"
#include "nrtk/atmos.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double pos[3];              /* VRS marker ECEF (m) */
    double elmin;               /* output elevation mask (rad) */
    int    sysmask;             /* output systems */
    int    pass_unfixed;        /* 1: output satellites without corrections with lock reset */
    double jump_thres;          /* correction jump (m) that triggers a lock reset */
    int    tide;                /* solid earth tide at the VRS */
    nrtk_sigplan_t plan[NRTK_NSYS];
} nrtk_vrs_opt_t;

void nrtk_vrs_opt_default(nrtk_vrs_opt_t *opt);

typedef struct {
    int    nsat;                /* output satellites */
    int    nlockreset;          /* lock time resets in this epoch */
    int    master;              /* master station */
    double dtr;                 /* master receiver clock used (s) */
} nrtk_vrs_info_t;

typedef struct nrtk_vrs nrtk_vrs_t;

nrtk_vrs_t *nrtk_vrs_create(const nrtk_vrs_opt_t *opt);
void nrtk_vrs_free(nrtk_vrs_t *v);
/* forget continuity (new VRS / master change): all lock times restart */
void nrtk_vrs_reset(nrtk_vrs_t *v);
const nrtk_vrs_opt_t *nrtk_vrs_getopt(const nrtk_vrs_t *v);

/* generate the VRS observations of one epoch. out: obsd_t[maxn], rcv=1. */
int nrtk_vrs_generate(nrtk_vrs_t *v, const nrtk_staepoch_t *master,
                      const nrtk_vcorr_t *corr, const nav_t *nav,
                      obsd_t *out, int maxn, nrtk_vrs_info_t *info);

/* degraded mode: pass through the observations of one station (prepared,
 * null antenna). Lock times restart when the station changes. */
int nrtk_vrs_passthrough(nrtk_vrs_t *v, const nrtk_staepoch_t *sta,
                         obsd_t *out, int maxn, nrtk_vrs_info_t *info);

#ifdef __cplusplus
}
#endif
#endif
