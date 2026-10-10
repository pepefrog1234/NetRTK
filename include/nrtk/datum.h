/*------------------------------------------------------------------------------
 * nrtk/datum.h : output datums - shift of the RTCM 1006 position from ITRF to
 *                a static legal datum (Taiwan TWD97[2010] / TWD97[2020])
 *
 * The network solution, the atmospheric interpolation and the virtual
 * observations stay in ITRF at the current epoch: the observations must match
 * the broadcast orbits, and a station coordinate error of 10 cm already shows
 * up as 15 cm in the double-difference residuals (0.8 L1 cycles). Taiwan moves
 * 3-7 cm/yr with a direction that changes across the island, so coordinates
 * frozen at 2010.0 / 2020.0 differ from ITRF(t) by 0.5-1.3 m / 0.2-0.5 m and
 * by 5-30 cm between the ends of one baseline. A legal datum can therefore
 * only enter where the rover sees it and the solver does not: the reference
 * position in RTCM 1006. The rover's result is "1006 position + baseline
 * vector"; shifting the 1006 position by delta moves the result by delta and
 * changes the baseline by |delta|*|b|/rho (0.3 mm at 5 km).
 *
 * delta is a vector field: for every station (and optional calibration point)
 * delta_i(t) = X_official_i - [X_ITRF_i(t0) + v_i (t - t0)], interpolated
 * linearly inside the Delaunay triangles of the points (barycentric weights,
 * continuous across triangle edges so a rover changing cells does not jump).
 * Outside the triangulation the plane of the nearest triangle is extended up
 * to extrap_max_km; farther away the nearest point's delta is used.
 * The field is evaluated once when a VRS is built (a session keeps its delta).
 *
 * Also here: the "theory chain" (velocity back to the datum epoch + IERS
 * ITRF2020 -> ITRFxx Helmert transformation) used to check the official
 * coordinates, and the TM2 projection used for display only.
 *----------------------------------------------------------------------------*/
#ifndef NRTK_DATUM_H
#define NRTK_DATUM_H
#include "nrtk/config.h"
#ifdef __cplusplus
extern "C" {
#endif

/* how the shift at a position was obtained */
enum { NRTK_DSHIFT_NONE=0,      /* no points: no shift available */
       NRTK_DSHIFT_TRI=1,       /* inside a triangle (barycentric) */
       NRTK_DSHIFT_EXTRAP=2,    /* outside: plane of the nearest triangle */
       NRTK_DSHIFT_NEAREST=3 }; /* nearest point's shift (far outside, or
                                   interp=nearest, or fewer than 3 points) */

typedef struct {
    int    mode;                /* NRTK_DSHIFT_xxx */
    int    pt[3];               /* point indices used (-1: unused) */
    double w[3];                /* weights */
    double extrap_km;           /* planar distance outside the triangulation (0 inside) */
    int    theory;              /* a point whose coordinates come from the theory chain was used */
    double enu[3];              /* delta in local east/north/up at the position (m) */
} nrtk_datum_info_t;

typedef struct nrtk_datum nrtk_datum_t;

/* build the shift field of configuration datum idatum from the stations that
 * have coordinates in it (enabled or not) and the calibration points.
 * Stations enabled in the network but without coordinates in this datum are
 * listed as "missing"; a missing station is "covered" when a calibration
 * point lies within cover_km of it. A network cell with an uncovered missing
 * station is datum-incomplete: the shift there could only be interpolated
 * across the station from farther points, so service in this datum is
 * refused in that cell (nrtk_datum_cell_incomplete). */
nrtk_datum_t *nrtk_datum_create(const nrtk_datum_def_t *def, int idatum,
                                const nrtk_station_t *sta, int nsta,
                                const nrtk_dpoint_t *pts, int npts,
                                int interp, double extrap_max_km, double cover_km);
nrtk_datum_t *nrtk_datum_from_config(const nrtk_config_t *cfg, int idatum);
void nrtk_datum_free(nrtk_datum_t *d);

const nrtk_datum_def_t *nrtk_datum_def(const nrtk_datum_t *d);
int  nrtk_datum_npoints(const nrtk_datum_t *d);
int  nrtk_datum_ntri(const nrtk_datum_t *d);
const char *nrtk_datum_point_id(const nrtk_datum_t *d, int i);
int  nrtk_datum_point_station(const nrtk_datum_t *d, int i);     /* station index, -1: calibration point */
const char *nrtk_datum_point_source(const nrtk_datum_t *d, int i);
/* ITRF position of point i at time t, and its shift delta_i(t) (ECEF, m) */
void nrtk_datum_point(const nrtk_datum_t *d, int i, gtime_t t, double *pos_itrf,
                      double *delta);
int  nrtk_datum_nmissing(const nrtk_datum_t *d);
int  nrtk_datum_missing(const nrtk_datum_t *d, int k);            /* station index */
int  nrtk_datum_missing_covered(const nrtk_datum_t *d, int k);    /* calibration point near it */
int  nrtk_datum_nuncovered(const nrtk_datum_t *d);                /* missing and not covered */
/* network cell (3 station indices, -1 unused) with an uncovered missing station */
int  nrtk_datum_cell_incomplete(const nrtk_datum_t *d, const int *tri);

/* shift at ITRF position pos (ECEF) at time t: X_datum = pos + delta.
 * returns the mode (NRTK_DSHIFT_NONE: delta is zero, nothing available) */
int  nrtk_datum_shift(const nrtk_datum_t *d, const double *pos, gtime_t t,
                      double *delta, nrtk_datum_info_t *info);
/* ITRF position of a position given in the datum (two iterations, < 1 mm) */
int  nrtk_datum_inverse(const nrtk_datum_t *d, const double *pos_datum, gtime_t t,
                        double *pos_itrf);
/* leave-one-out: shift at point i interpolated from the other points.
 * returns the mode (0 if fewer than 3 other points) */
int  nrtk_datum_loo(const nrtk_datum_t *d, int i, gtime_t t, double *delta,
                    nrtk_datum_info_t *info);

/* theory chain --------------------------------------------------------------*/
/* IERS 14-parameter transformation ITRF2020 -> frame at epoch (decimal year):
 * X_frame = X + T + D X + R x X, parameters of the IERS ITRF product centre
 * (Transfo-ITRF2020_TRFs.txt, reference epoch 2015.0). frame: ITRF2014,
 * ITRF2008, ITRF2005, ITRF2000, ITRF97, ITRF96, ITRF94. returns 0 if unknown */
int  nrtk_helmert_itrf2020(const char *frame, double epoch, const double *x,
                           double *y);
/* position (ITRF2020, at decimal year t0, velocity m/yr) moved back to the
 * datum epoch and transformed to the datum's ITRF frame; coseismic steps
 * between the epochs are not modelled. returns 0 if the frame is unknown */
int  nrtk_datum_theory(const double *pos, double t0, const double *vel,
                       const char *frame, double epoch, double *out);

/* TM2 projection (display only) ----------------------------------------------
 * Transverse Mercator, GRS80, k0 = 0.9999, false easting 250 000 m, central
 * meridian lon0_deg: 121 (Taiwan, EPSG:3826) or 119 (Penghu, Kinmen, Matsu,
 * EPSG:3825). pos: lat, lon (rad). Krueger series to n^4: < 1 mm in Taiwan. */
void nrtk_tm2(const double *pos, double lon0_deg, double *E, double *N);
double nrtk_tm2_zone(double lon_rad);     /* 121 or 119 */

#ifdef __cplusplus
}
#endif
#endif
