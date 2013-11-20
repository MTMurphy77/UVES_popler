/******************************************************************************
AST_VORBIT: Radial velocity component of the Earth-Moon barycenter relative
to the Sun.

double	ra		# Right ascension of observation (hours)
double	dec		# Declination of observation (degrees)
double	epoch		# Julian epoch of observation
double	v		# Component of orbital velocity (km/s)

******************************************************************************/

#include <math.h>
#include "astron.h"
#include "error.h"

double ast_vorbit(double ra, double dec, double epoch) {

  double t=0.0,manom=0.0,lperi=0.0,oblq=0.0,eccen=0.0,tanom=0.0,slong=0.0;
  double r=0.0,d=0.0,l=0.0,b=0.0,vorb=0.0;
  double eccensq=0.0,eccencu=0.0;
  
  /* T is the number of Julian centuries since J1900 */
  t=(ast_epoch2jd(epoch)-2415020.0)/36525.0;
  /* MANOM is the mean anomaly of the Earth's orbit (degrees)
     LPERI is the mean longitude of perihelion (degrees)
     OBLQ is the mean obliquity of the ecliptic (degrees)
     ECCEN is the eccentricity of the Earth's orbit (dimensionless) */
  manom=358.47583+t*(35999.04975-t*(0.000150+t*0.000003));
  lperi=101.22083+t*(1.7191733+t*(0.000453+t*0.000003));
  oblq=23.452294-t*(0.0130125+t*(0.00000164-t*0.000000503));
  eccen=0.01675104-t*(0.00004180+t*0.000000126);
  eccensq=eccen*eccen; eccencu=eccensq*eccen;

  /* Convert to principle angles */
  manom=AST_MOD(manom,360.0); lperi=AST_MOD(lperi,360.0);
  /* Convert to radians */
  r=ra*15.0*C_RPDEG; d=dec*C_RPDEG;
  manom*=C_RPDEG; lperi*=C_RPDEG; oblq*=C_RPDEG;

  /* TANOM is the true anomaly (approximate formula) (radians) */
  tanom=manom+(2.0*eccen-0.25*eccencu)*sin(manom)+1.25*eccensq*
    sin(2.0*manom)+13.0/12.0*eccencu*sin(3.0*manom);

  /* SLONG is the true longitude of the Sun seen from the Earth (radians) */
  slong=lperi+tanom+C_PI;

  /* L and B are the longitude and latitude of the star in the orbital
     plane of the Earth (radians) */
  if (!ast_coord(0.0,0.0,-C_PI_2,C_PI_2-oblq,r,d,&l,&b))
    errormsg("ast_vorbit(): Unknown error returned from ast_coord()");

  /* VORB is the component of the Earth's orbital velocity perpendicular to
     the radius vector (km/s) where the Earth's semi-major axis is
     149598500 km and the year is 365.2564 days */
  vorb=((C_2PI/365.2564)*149598500.0/sqrt(1.0-eccensq))/86400.0;

  /* V is the projection onto the line of sight to the observation of the
     velocity of the Earth-Moon barycenter with respect to the Sun
     (km/s) */
  vorb*=cos(b)*(sin(slong-l)-eccen*sin(lperi-l));

  return vorb;

}
