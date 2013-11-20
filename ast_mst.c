/******************************************************************************
AST_MST: Mean sidereal time of the epoch at the given longitude. This
procedure may be used to obtain Greenwich Mean Sidereal Time (GMST) by
setting the longitude to 0.

double	epoch		# Epoch
double	lon     	# Longitude in degrees

******************************************************************************/

#include "astron.h"

double ast_mst(double epoch, double lon) {

  double jd=0.0,ut=0.0,t=0.0,st=0.0;

  /* Determine JD and UT, and T (JD in centuries from J2000.0) */
  jd=ast_epoch2jd(epoch); ut=(jd-(int)jd-0.5)*24.0;
  t=(jd-C_JD2000)/(100.0*C_JYEAR);

  /* The GMST at 0 UT in seconds is a power series in T */
  st=24110.54841+t*(8640184.812866+t*(0.093104-t*6.2e-6));

  /* Correct for longitude and convert to standard hours */
  st=AST_MOD((st/C_SECHR+ut-lon/15.0),24.0);
  if (st < 0.0) st+=24.0;

  return st;

}
