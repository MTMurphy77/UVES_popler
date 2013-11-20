/****************************************************************************
* Definition of the Edlen formula for finding the vacuum wavelength from a
* given air wavelength. The formula for taking a given vacuum wavelength
* and converting it into an air wavelength is given in Edlen B., 1966,
* Metrologia, 2, 71 and is the formula used by Palmer B. & Engleman R. Jr.,
* 1983, Atlas of the thorium spectrum, Los Alamos National Laboratory, Los
* Alamnos, NM and Norlen G., 1973, Phys. Scr., 8, 249 to convert the vacuum
* wavelengths of the ThAr spectrum to air wavelengths. Therefore, one must
* solve the formula below iteratively to arrive at the correct vacuum
* wavelength. If one does not iterate, the errors in the refractive index
* will be in error by as much as 4x10^{-9} at 10,000 Angstroems (smaller at
* shorter optical wavelengths). See Murphy M. et al., 2001, MNRAS, 327,
* 1223 for more discussion.
*
* The routine below solves the Edlen formula iteratively, returning the
* correct vacuum wavelength to the precision specified by EDLEN_A2V_TOL (in
* Angstroems). The argument, airwl, must be in Angstroems.
****************************************************************************/

#include <math.h>
#include "error.h"

#define EDLEN_A2V_TOL   1.e-14 /* Precision of refractive index conversion */
#define EDLEN_A2V_ITER 20      /* Maximum number of iterations allowed     */

double edlen_a2v(double airwl) {

  double n=1.0,n_prev=0.0; /* First guess refractive index is 1.0 */
  double sigmasq=0.0;
  int    iter=0;

  if (airwl<=0.0) {
    nferrormsg("edlen_a2v(): Input wavelength %lf invalid",airwl); 
    return 0.0;
  }

  while (fabs(n-n_prev)>EDLEN_A2V_TOL && iter<EDLEN_A2V_ITER) {
    n_prev=n; sigmasq=1.e8/(airwl*airwl*n_prev*n_prev);
    n=1.e-8*(15997.0/(38.90-sigmasq)+2406030.0/(130.0-sigmasq)+8342.13)+1.0;
    iter++;
  }
  if (iter>=EDLEN_A2V_ITER) {
    nferrormsg("edlen_a2v(): Did not converge on vacuum wavelength\n\
\tin %d iterations to refractive index precision %e\n\
\tfor air wavelength %lf",EDLEN_A2V_ITER,EDLEN_A2V_TOL,airwl);
    return 0.0;
  }

  return airwl*n;

}
