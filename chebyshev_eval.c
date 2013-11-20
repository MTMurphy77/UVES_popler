/*************************************************************************** 
CHEBYSHEV_EVAL: Evalutate the Chebyshev polynomial described by the
coefficients array p at the value of x provided. THE VALUE OF x INPUT
MUST BE A RENORMALIZATION OF ONE'S OWN WORKING VALUES: x =
[2.0*x_own-b-a]/[b-a].
****************************************************************************/

#include "error.h"

double chebyshev_eval(double x, double *coef, int ord) {

  double twox=0.0,d=0.0,dd=0.0,sv=0.0;
  int    i=0;

  if (x<-1.0 || x>1.0)
    errormsg("chebyshev_eval(): Input value of x (=%lg)\n\
\tmust be >=-1.0 and <=1.0",x);

  twox=2.0*x;
  for (i=ord-1; i>0; i--) { sv=d; d=twox*d-dd+coef[i]; dd=sv; }

  return x*d-dd+coef[0];
  
}
