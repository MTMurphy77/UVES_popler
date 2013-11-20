/*************************************************************************** 
SVDFIT_LEGENDRE: This is the FLEG algorithm taken from Numercial Recipes.

I have modified the routine in the following ways:

0) Double precision

****************************************************************************/

int svdfit_legendre(double x, double *pl, int nl) {

  int    i=0;
  double twox=0.0,f2=0.0,f1=0.0,d=1.0;
  
  pl[0]=1.0;
  if (nl>1) pl[1]=x;
  if (nl>2) {
    twox=2.0*x; f2=x;
    for (i=2; i<nl; i++) {
      f1=d++; f2+=twox;
      pl[i]=(f2*pl[i-1]-f1*pl[i-2])/d;
    }
  }

  return 1;

}
