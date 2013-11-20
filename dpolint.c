/***************************************************************************** 
POLINT: This is the POLINT algorithm taken from Numerical Recipes. Here's
what NR has to say about it:

"Given arrays xa[1..n] and ya[1..n], and given a value x, this routine
returns a value y, and an error estimate dy. If P(x) is the polynomial
of degree N-1 such that P(x_a_i) = y_a_i, i = 1, . . . , n, then the
returned value y = P(x)."

I have made the following changes:
0) Double precision inputs and outputs and returns an interger
1) Own error handling routines.
*****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include "memory.h"
#include "error.h"

int dpolint(double *xa, double *ya, int n, double x, double *y, double *dy) {

  int       i=0,m=0,ns=0;
  double    den=0.0,dif=0.0,dift=0.0,ho=0.0,hp=0.0,w=0.0;
  double    *c=NULL,*d=NULL;

  /* Allocate temporary memory */
  if ((c=darray(n))==NULL) {
    nferrormsg("dpolint(): Cannot allocate memory to c array of size %d",n);
    return 0;
  }
  if ((d=darray(n))==NULL) {
    nferrormsg("dpolint(): Cannot allocate memory to d array of size %d",n);
    return 0;
  }

  dif=fabs(x-xa[0]);
  for (i=0; i<n; i++) {
    if ((dift=fabs(x-xa[i]))<dif) { ns=i; dif=dift; }
    c[i]=ya[i]; d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=0; m<n-1; m++) {
    for (i=0; i<n-m-1; i++) {
      ho=xa[i]-x; hp=xa[i+m+1]-x; w=c[i+1]-d[i];
      if ((den=ho-hp)==0.0) {
        nferrormsg("dpolint(): Attempt to divide by zero.\n\
\tInput x array may contain two or more identical entries"); return 0;
      }
      den=w/den; d[i]=hp*den; c[i]=ho*den;
    }
    *y+=(*dy=(2*(ns+1)<(n-m-1) ? c[ns+1] : d[ns--]));
  }

  /* Clean up */
  free(d); free(c);

  return 1;

}
