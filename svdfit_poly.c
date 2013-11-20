/*************************************************************************** 
SVDFIT_POLY: This is the FPOLY algorithm taken from Numercial Recipes.

I have modified the routine in the following ways:

0) Double precision

****************************************************************************/

int svdfit_poly(double x, double *p, int np) {

  int i=0;

  p[0]=1.0;
  for (i=1; i<np; i++) p[i]=p[i-1]*x;

  return 1;

}
