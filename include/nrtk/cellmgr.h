/*------------------------------------------------------------------------------
 * nrtk/cellmgr.h : station network cells (Delaunay triangulation), master /
 *                  auxiliary station selection, station health
 *----------------------------------------------------------------------------*/
#ifndef NRTK_CELLMGR_H
#define NRTK_CELLMGR_H
#include "nrtk/common.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int    valid;               /* a cell could be formed */
    int    master;              /* master station index */
    int    naux;
    int    aux[NRTK_MAXSTA];    /* auxiliary stations (triangle partners first) */
    int    tri[3];              /* containing triangle, -1 if none */
    int    extrap;              /* V outside the triangulation (extrapolation) */
    int    reject;              /* V farther than max distance from any station */
    int    nearest;             /* nearest healthy station */
    double dist_km;             /* distance V - master (km) */
} nrtk_cell_t;

typedef struct {
    double max_dist_km;         /* reject service beyond this distance (70) */
    double lsm_radius_km;       /* extra aux stations within this radius (50) */
    int    max_aux;             /* max aux stations per cell */
    double max_baseline_km;     /* compute baselines up to this length */
} nrtk_cellmgr_opt_t;

void nrtk_cellmgr_opt_default(nrtk_cellmgr_opt_t *opt);

typedef struct nrtk_cellmgr nrtk_cellmgr_t;

nrtk_cellmgr_t *nrtk_cellmgr_create(const nrtk_cellmgr_opt_t *opt,
                                    const nrtk_station_t *sta, int nsta);
void nrtk_cellmgr_free(nrtk_cellmgr_t *cm);

/* station health; changing health re-triangulates the healthy stations.
 * returns 1 if the triangulation changed */
int  nrtk_cellmgr_set_health(nrtk_cellmgr_t *cm, int ista, int healthy);
int  nrtk_cellmgr_healthy(const nrtk_cellmgr_t *cm, int ista);
/* select the cell for virtual position vpos (ECEF). returns cell.valid */
int  nrtk_cellmgr_select(const nrtk_cellmgr_t *cm, const double *vpos,
                         nrtk_cell_t *cell);
/* current triangles of healthy stations */
int  nrtk_cellmgr_triangles(const nrtk_cellmgr_t *cm, int (*tri)[3], int max);
/* station pairs (i<j) for which baselines are computed */
int  nrtk_cellmgr_baselines(const nrtk_cellmgr_t *cm, int (*pair)[2], int max);
/* nearest (healthy) station to vpos, -1 if none */
int  nrtk_cellmgr_nearest(const nrtk_cellmgr_t *cm, const double *vpos,
                          int healthy_only, double *dist_km);

/* planar Delaunay triangulation (exposed for tests). x,y: n points.
 * returns number of triangles written to tri (counter-clockwise). */
int  nrtk_delaunay(const double *x, const double *y, int n, int (*tri)[3],
                   int max);

#ifdef __cplusplus
}
#endif
#endif
