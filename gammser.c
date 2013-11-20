/****************************************************************************
GAMMSER: This is the incomplete gamma function algorithm taken from Numerical
Recipes. Here's what NR has to say;

"Returns the incomplete gamma function P(a,x) evaluated by its series
representation as gser. Also returns lnGamma(a) as glnx".

Best to use the higher level function gammq or gammp to calculate this.

I have changed the raw routine to include the following:
0) Double precision is used
1) Returns an integer

****************************************************************************/

#include <math.h>
#include "gamm.h"
#include "error.h"

int gammser(double *gser, double a, double x, double *gln) {

  double sum,del,ap;
  int    n;
  
  *gln=gammln(a);
  if (x<=0.0) {
    if (x<0.0) {
      warnmsg("gammser(): x (=%lf) is less than 0, which is invalid.",x);
      return 0;
    }
    *gser=0.0;
    return 1;
  }
  else {
    ap=a;
    del=sum=1.0/a;
    for (n=1; n<GAMM_ITMAX; n++) {
      ++ap;
      del*=x/ap;
      sum+=del;
      if (fabs(del)<fabs(sum)*GAMM_EPS) {
	*gser=sum*exp(-x+a*log(x)-(*gln));
	return 1;
      }
    }
    warnmsg("gammser(): a (=%lf) too large, GAMM_ITMAX (=%d) too small.\n\
\tTry increasing ITMAX in gamm.h.",a,GAMM_ITMAX);
    return 0;
  }
}
