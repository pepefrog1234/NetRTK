/*------------------------------------------------------------------------------
 * nrtk/tropo.h : troposphere model (Saastamoinen zenith delays + GMF)
 *----------------------------------------------------------------------------*/
#ifndef NRTK_TROPO_H
#define NRTK_TROPO_H
#include "nrtk/common.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Global Mapping Function (Boehm 2006). hgt: height above MSL (m) */
int nrtk_gmf(double mjd, double lat, double lon, double hgt, double zd,
             double *gmfh, double *gmfw);

/* Saastamoinen zenith hydrostatic/wet delays with a standard atmosphere
 * (relative humidity humi). pos: geodetic {lat,lon,h} (rad,rad,m) */
void nrtk_trop_zenith(const double *pos, double humi, double *zhd, double *zwd);

/* slant tropospheric delay model (m): Saastamoinen zenith delays mapped with
 * GMF. Used both for the network residual reduction and at the VRS. */
double nrtk_trop_slant(gtime_t t, const double *pos, const double *azel);

#ifdef __cplusplus
}
#endif
#endif
