/*------------------------------------------------------------------------------
 * nrtk/antenna.h : receiver antenna phase centre model (ANTEX 1.4)
 *
 * Unlike RTKLIB's reader (GPS frequencies only, elevation-only PCV), this
 * keeps every system/frequency block (G01, G02, G05, E01, E05, E07, E08, E06,
 * C01, C02, C05, C06, C07, J01, ...) and azimuth-dependent PCV.
 *----------------------------------------------------------------------------*/
#ifndef NRTK_ANTENNA_H
#define NRTK_ANTENNA_H
#include "nrtk/common.h"
#ifdef __cplusplus
extern "C" {
#endif

#define NRTK_ANT_MAXFRQ 48

typedef struct {
    char    code[4];          /* ANTEX frequency code, e.g. "G01" */
    double  neu[3];           /* PCO north/east/up (m) */
    double  zen1,zen2,dzen;   /* zenith grid (deg) */
    double  dazi;             /* azimuth step (deg), 0: no azimuth dependence */
    int     nzen,nazi;
    double *noazi;            /* [nzen] PCV (m) */
    double *azi;              /* [nazi][nzen] PCV (m), azimuth 0..360 */
} nrtk_antfrq_t;

typedef struct {
    char          type[24];   /* antenna type (20 chars, trimmed) */
    char          radome[8];  /* radome (4 chars, trimmed) */
    int           nf;
    nrtk_antfrq_t f[NRTK_ANT_MAXFRQ];
} nrtk_ant_t;

typedef struct {
    int n,nmax;
    nrtk_ant_t *ant;
} nrtk_atx_t;

int  nrtk_atx_read(const char *file, nrtk_atx_t *atx);   /* receiver antennas */
void nrtk_atx_free(nrtk_atx_t *atx);
/* find antenna by "TYPE RADOME" (radome optional; falls back to NONE) */
const nrtk_ant_t *nrtk_atx_find(const nrtk_atx_t *atx, const char *name);

/* range correction (m) for system sys, carrier freq (Hz), azel at station:
 * modelled observation = range(marker) + dr. Includes ARP eccentricity del
 * (E/N/U). ant may be NULL (only eccentricity is applied). */
double nrtk_ant_dr(const nrtk_ant_t *ant, int sys, double freq,
                   const double *azel, const double *del, int pcv);

#ifdef __cplusplus
}
#endif
#endif
