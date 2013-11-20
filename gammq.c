/****************************************************************************
GAMMQ: This is the incomplete gamma function algorithm taken from Numerical
Recipes. Here's what NR has to say;

"Returns the incomplete gamma function Q(a,x)=1-P(a,x)".

Remember that if you are trying to calculate the probability Q that
chi-squared should exceed a particular value of chi-squared by chance then
you want Q = gammq(0.5*NDF, 0.5*chi-squared).

I have changed the raw routine to include the following:
0) Double precision is used

****************************************************************************/

#include "gamm.h"
#include "error.h"

double gammq(double a, double x) {

  double gser,gcf,gln;

  if (x<0.0) errormsg("gammq(): Invalid x-value.");
  if (a<=0.0) errormsg("gammq(): Invalid a-value.");
  if (x<(a+1.0)) {
    /* Use series representation */
    if (!gammser(&gser,a,x,&gln))
      errormsg("gammq(): Error returned from gammser().");
    return 1.0-gser;
  }
  else {
    /* Use the continued fraction representation */
    if (!gammcf(&gcf,a,x,&gln))
      errormsg("gammq() Error returned from gammcf().");
    return gcf;
  }
}

