/*------------------------------------------------------------------------------
 * tropo.c : Saastamoinen zenith delays (standard atmosphere) mapped with GMF
 *----------------------------------------------------------------------------*/
#include <math.h>
#include "nrtk/tropo.h"

extern void nrtk_trop_zenith(const double *pos, double humi, double *zhd,
                             double *zwd)
{
    double hgt,pres,temp,e;

    hgt=pos[2]-geoidh(pos);             /* height above mean sea level */
    if (hgt<0.0) hgt=0.0;
    if (hgt>1E4) hgt=1E4;
    pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
    temp=15.0-6.5E-3*hgt+273.16;
    e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));
    *zhd=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3);
    *zwd=0.002277*(1255.0/temp+0.05)*e;
}
extern double nrtk_trop_slant(gtime_t t, const double *pos, const double *azel)
{
    const double ep2000[]={2000,1,1,12,0,0};
    double zhd,zwd,mjd,mh,mw;

    if (azel[1]<=0.0||pos[2]<-1000.0||pos[2]>2E4) return 0.0;
    nrtk_trop_zenith(pos,0.7,&zhd,&zwd);
    mjd=51544.5+timediff(t,epoch2time(ep2000))/86400.0;
    nrtk_gmf(mjd,pos[0],pos[1],pos[2]-geoidh(pos),PI/2.0-azel[1],&mh,&mw);
    return zhd*mh+zwd*mw;
}
