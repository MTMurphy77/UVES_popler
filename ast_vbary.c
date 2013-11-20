/******************************************************************************
AST_VBARY: Radial velocity component of center of the Earth relative to the
barycenter of the Earth-Moon system.

double	ra		# Right ascension of observation (hours)
double	dec		# Declination of observation (degrees)
double	epoch		# Julian epoch of observation
double	v		# Component of orbital velocity (km/s)

******************************************************************************/

#include <math.h>
#include "astron.h"
#include "error.h"

double ast_vbary(double ra, double dec, double epoch) {

  double t=0.0,oblq=0.0,omega=0.0,llong=0.0,lperi=0.0,inclin=0.0,em=0.0;
  double anom=0.0,vmoon=0.0;
  double r=0.0,d=0.0,l=0.0,b=0.0,lm=0.0,bm=0.0;
  double emsq=0.0,emcu=0.0;

  /* T is the number of Julian centuries since J1900 */
  t=(ast_epoch2jd(epoch)-2415020.0)/36525.0;

  /*
    OBLQ is the mean obliquity of the ecliptic
    OMEGA is the longitude of the mean ascending node
    LLONG is the mean lunar longitude (should be 13.1763965268)
    LPERI is the mean lunar longitude of perigee
    INCLIN is the inclination of the lunar orbit to the ecliptic
    EM is the eccentricity of the lunar orbit (dimensionless)
    All quantities except the eccentricity are in degrees */

  oblq=23.452294-t*(0.0130125+t*(0.00000164-t*0.000000503));
  omega=259.183275-t*(1934.142008+t*(0.002078+t*0.000002));
  llong=270.434164+t*(481267.88315+t*(-0.001133+t*0.0000019))-omega;
  lperi=334.329556+t*(4069.034029-t*(0.010325+t*0.000012))-omega;
  em=0.054900489; inclin=5.1453964;
  emsq=em*em; emcu=emsq*em;

  /* Determine true longitude.  Compute mean anomaly, convert to true
     anomaly (approximate formula), and convert back to longitude. The mean
     anomaly is only approximate because LPERI should be the true rather
     than the mean longitude of lunar perigee */
  lperi*=C_RPDEG; llong*=C_RPDEG; anom=llong-lperi;
  anom+=(2.0*em-0.25*emcu)*sin(anom)+1.25*emsq*sin(2.0*anom)+
    13.0/12.0*emcu*sin(3.0*anom);
  llong=anom+lperi;

  /* L and B are the ecliptic longitude and latitude of the observation. LM
     and BM are the lunar longitude and latitude of the observation in the
     lunar orbital plane relative to the ascending node */
  r=C_RPDEG*ra*15.0; d=C_RPDEG*dec; omega*=C_RPDEG; oblq*=C_RPDEG;
  inclin*=C_RPDEG;

  if (!ast_coord(0.0,0.0,-C_PI_2,C_PI_2-oblq,r,d,&l,&b))
    errormsg("ast_vbary(): Unknown error returned from ast_coord()");
  if (!ast_coord(omega,0.0,omega-C_PI_2,C_PI_2-inclin,l,b,&lm,&bm))
    errormsg("ast_vbary(): Unknown error returned from ast_coord()");

  /* VMOON is the component of the lunar velocity perpendicular to the
     radius vector.  V is the projection onto the line of sight to the
     observation of the velocity of the Earth's center with respect to the
     Earth-Moon barycenter.  The 81.53 is the ratio of the Earth's mass to
     the Moon's mass */
  vmoon=(C_2PI/27.321661)*384403.12040/sqrt(1.0-emsq)/86400.0;
  vmoon*=cos(bm)*(sin(llong-lm)-em*sin(lperi-lm))/81.53;

  return vmoon;

}
