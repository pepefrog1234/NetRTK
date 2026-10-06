/*------------------------------------------------------------------------------
 * nrtk/ingest.h : data ingest - RINEX files (offline) and RTCM 3 streams
 *                 (real time), decoded with RTKLIB into obsd_t / nav_t
 *----------------------------------------------------------------------------*/
#ifndef NRTK_INGEST_H
#define NRTK_INGEST_H
#include <stdio.h>
#include "nrtk/common.h"
#ifdef __cplusplus
extern "C" {
#endif

/* offline: epoch-wise RINEX observation reader (constant memory) ------------*/
typedef struct {
    FILE    *fp;
    rnxctr_t rnx;
    int      rcv;
    int      buffered;          /* an epoch is buffered in rnx.obs */
    int      eof;
    char     path[1024];
    char     tmp[1024];         /* uncompressed temporary file, if any */
} nrtk_rnxsrc_t;

int     nrtk_rnxsrc_open (nrtk_rnxsrc_t *s, const char *file, int rcv);
void    nrtk_rnxsrc_close(nrtk_rnxsrc_t *s);
/* time of the next buffered epoch ({0} at end of file) */
gtime_t nrtk_rnxsrc_peek (nrtk_rnxsrc_t *s);
/* observations at time t (|dt|<DTTOL); older epochs are skipped.
 * returns number of observations, 0 if none at t, -1 at end of file */
int     nrtk_rnxsrc_read (nrtk_rnxsrc_t *s, gtime_t t, obsd_t *obs, int max);
/* RINEX header station info (after open) */
const sta_t *nrtk_rnxsrc_sta(const nrtk_rnxsrc_t *s);

/* navigation files: blank/comma separated list, RTKLIB wildcards allowed */
int     nrtk_nav_read(const char *files, nav_t *nav);
void    nrtk_nav_free(nav_t *nav);
/* deep copy of broadcast ephemerides (snapshot for one processing epoch) */
int     nrtk_nav_copy(nav_t *dst, const nav_t *src);
/* merge one ephemeris into nav (newest per satellite and IODE kept) */
int     nrtk_nav_add_eph(nav_t *nav, const eph_t *eph);

/* real time: RTCM 3 stream input ------------------------------------------*/
typedef struct {
    stream_t str;
    rtcm_t   rtcm;
    int      opened;
    char     url[512];
    char     msg[256];
    long     nbyte;             /* received bytes */
    long     nepoch;            /* decoded observation epochs */
    gtime_t  last_rx;           /* wall clock of last epoch (GPST) */
} nrtk_strsrc_t;

/* url: ntrip://user:pass@host:port/mount, tcpcli://host:port,
 *      tcpsvr://:port, file:///path[::T[::x<speed>]], serial://dev:baud */
int  nrtk_strsrc_open (nrtk_strsrc_t *s, const char *url);
void nrtk_strsrc_close(nrtk_strsrc_t *s);
/* read and decode available bytes. For each complete observation epoch
 * cb_obs is called; for each ephemeris cb_eph. returns bytes read, -1 error */
int  nrtk_strsrc_poll (nrtk_strsrc_t *s,
                       void (*cb_obs)(void *ctx, const obsd_t *obs, int n),
                       void (*cb_eph)(void *ctx, const eph_t *eph),
                       void (*cb_sta)(void *ctx, const sta_t *sta),
                       void *ctx);
/* convert URL to RTKLIB stream type/path. returns STR_xxx or 0 */
int  nrtk_url2str(const char *url, char *path, int size);

#ifdef __cplusplus
}
#endif
#endif
