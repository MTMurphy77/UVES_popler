/*************************************************************************** 
SVDFIT_CHEBYSHEV: Calculates the Chebyshev series for use with the
svdfit() fitting algorithm. THE VALUE OF x INPUT MUST BE A
RENORMALIZATION OF ONE'S OWN WORKING VALUES: x = [2.0*x_own-b-a]/[b-a].
****************************************************************************/

#include "error.h"

int svdfit_chebyshev(double x, double *pt, int nt) {

  int    i=0;
  double twox=0.0;

  if (x<-1.0 || x>1.0)
    errormsg("svdfit_chebyshev(): Input value of x (=%lf)\n\
\tmust be >=-1.0 and <=1.0",x);

  pt[0]=1.0;
  if (nt>1) pt[1]=x;
  if (nt>2) {
    twox=2.0*x;
    for (i=2; i<nt; i++) pt[i]=twox*pt[i-1]-pt[i-2];
  }

  return 1;

}
