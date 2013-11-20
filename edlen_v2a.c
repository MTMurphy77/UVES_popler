/****************************************************************************
* Definition of the Edlen formula. The formula for taking a given vacuum
* wavelength and converting it into an air wavelength is given in Edlen B.,
* 1966, Metrologia, 2, 71 and is the formula used by Palmer B. & Engleman
* R. Jr., 1983, Atlas of the thorium spectrum, Los Alamos National
* Laboratory, Los Alamnos, NM and Norlen G., 1973, Phys. Scr., 8, 249 to
* convert the vacuum wavelengths of the ThAr spectrum to air
* wavelengths. The argument, airwl, must be in Angstroems.
****************************************************************************/

#include "error.h"

double edlen_v2a(double vacwl) {

  double n=0.0;
  double sigmasq=0.0;

  if (vacwl<=0.0) {
    nferrormsg("edlen_v2a(): Input wavelength %lf invalid",vacwl); 
    return 0.0;
  }

  sigmasq=1.e8/(vacwl*vacwl);
  n=1.e-8*(15997.0/(38.90-sigmasq)+2406030.0/(130.0-sigmasq)+8342.13)+1.0;

  return vacwl/n;

}
