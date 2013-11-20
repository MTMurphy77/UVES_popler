/*************************************************************************** 
LEGENDRE_EVAL: Evalutate the Legendre polynomial described by the
coefficients array p at the value of x provided.
****************************************************************************/

#include "fit.h"

double legendre_eval(double x, double *coef, int ord) {

  double y=0.0;
  double pl[ord];
  int i=0;

  /* Find calculate the values of the Legendre polynomials at the
     input value of x */
  svdfit_legendre(x,pl,ord);

  /* Evaluate the Legendre expansion series */
  for (i=ord-1; i>=0; i--) y+=coef[i]*pl[i];

  return y;

}
