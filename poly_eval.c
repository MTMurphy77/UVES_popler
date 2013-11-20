/*************************************************************************** 
POLY_EVAL: Evalutate the polynomial described by the coefficients
array p at the value of x provided.
****************************************************************************/

double poly_eval(double x, double *coef, int ord) {

  double y=0.0,power=1.0;
  int i=0;

  for (i=0; i<ord; i++) {
    if (i) power*=x; y+=coef[i]*power;
  }

  return y;

}
