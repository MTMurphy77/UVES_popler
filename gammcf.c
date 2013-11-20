/****************************************************************************
GAMMCF: This is the incomplete gamma function algorithm taken from Numerical
Recipes. Here's what NR has to say;

"Returns the incomplete gamma function Q(a,x) evaluated by its continued
fraction representation as gcf. Also returns lnGamma(a) as glnx".

Best to use the higher level function gammq or gammp to calculate this.

I have changed the raw routine to include the following:
0) Double precision is used
1) Returns an integer

****************************************************************************/

#include <math.h>
#include "gamm.h"
#include "error.h"

int gammcf(double *gcf, double a, double x, double *gln) {

  double an,b,c,d,del,h;
  int    i;
  
  *gln=gammln(a);
  b=x+1.0-a; c=1.0/GAMM_FPMIN; d=1.0/b; h=d;

  for (i=1; i<=GAMM_ITMAX; i++) {
    an=-(double)i*((double)i-a); b+=2.0; d=an*d+b;
    if (fabs(d)<GAMM_FPMIN) d=GAMM_FPMIN;
    c=b+an/c;
    if (fabs(c)<GAMM_FPMIN) c=GAMM_FPMIN;
    d=1.0/d; del=d*c; h*=del;
    if (fabs(del-1.0)<GAMM_EPS) break;
  }
  if (i>GAMM_ITMAX) {
    warnmsg("gammcf(): a (=%lf) too large, GAMM_ITMAX (=%d) too small.\n\
\tTry increasing GAMM_ITMAX in gamm.h.",a,GAMM_ITMAX);
    return 0;
  }

  *gcf=exp(-x+a*log(x)-(*gln))*h;
  
  return 1;
}
