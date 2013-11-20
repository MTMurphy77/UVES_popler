/******************************************************************************
AST_VROTATE: Radial velocity component of the observer relative to the
center of the Earth due to the Earth's rotation.

double	ra		# Right Ascension of observation (hours)
double	dec		# Declination of observation (degrees)
double	epoch		# Epoch of observation (Julian epoch)
double	lat	        # Latitude (degrees)
double	lon       	# Longitude (degrees; West is positive)
double	alt     	# Altitude (meters)
double	v		# Velocity (km / s)

******************************************************************************/

#include <math.h>
#include "astron.h"

double ast_vrotate(double ra, double dec, double epoch, double lat, double lon,
		   double alt) {

  double  latd=0.0,dlatd=0.0,r=0.0,vc=0.0,lmst=0.0;
  
  /* LATD is the latitude in radians */
  latd = C_RPDEG*lat;

  /* Reduction of geodetic latitude to geocentric latitude. Dlatd is in
     arcseconds */
  dlatd=-(11.0*60.0+32.743000)*sin(2.0*latd)+1.163300*sin(4.0*latd)-
    0.002600*sin(6.0*latd);
  latd+=C_RPDEG*dlatd/3600.0;

  /* R is the radius vector from the Earth's center to the observer
     (meters).  Vc is the corresponding circular velocity (meters/sidereal
     day converted to km / sec). (sidereal day = 23.934469591229 hours
     (1986)) */
  r=6378160.0*(0.998327073+0.00167643800*cos(2.0*latd)-0.00000351*cos(4.0*latd)
	       +0.000000008*cos(6.0*latd))+alt;
  vc=C_2PI*(r/1000.0)/(23.934469591229*3600.0);

  /* Project the velocity onto the line of sight to the star */
  lmst=ast_mst(epoch,lon);
  vc*=cos(latd)*cos(C_RPDEG*dec)*sin(C_RPDEG*(ra-lmst)*15.0);

  return vc;

}
