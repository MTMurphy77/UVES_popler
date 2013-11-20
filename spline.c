/*************************************************************************** 
SPLINE: This is the spline algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi
= f(xi), with x1 < x2 < .. . < xN, and given values yp1 and ypn for the
first derivative of the interpolating function at points 1 and n,
respectively, this routine returns an array y2[1..n] that contains the
second derivatives of the interpolating function at the tabulated points
xi. If yp1 and/or ypn are equal to 1 × 10^30 or larger, the routine is
signaled to set the corresponding boundary condition for a natural spline,
with zero second derivative on that boundary."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.
2) Arrays are 0-indexed
3) Error checking done with own error codes
4) Option for doing natural spline now encoded into "opt" argument. If
   opt=0 do normal spline and use given values of derivatives, if opt=1 do
   natural spline.
****************************************************************************/

#include <stdlib.h>
#include "memory.h"
#include "error.h"

int spline(double *x, double *y, int n, double yp1, double ypn, double *y2,
	   int opt) {

  double p=0.0,qn=0.0,sig=0.0,un=0.0;
  double *u=NULL;
  int    i=0,k=0;

  /* Check to make sure array length is reasonable */
  if (n<2) errormsg("spline(): Number of data points (=%d) must\n\
\tbe greater than 1",n);

  /* Allocate memory for temporary array */
  if ((u=darray(n-1))==NULL)
    errormsg("spline(): Could not allocate memory for u\n\
\tarray of size %d",n-1);

  if (opt) y2[0]=u[0]=0.0;
  else { y2[0]=-0.5; u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1); }
  for (i=1; i<n-1; i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]); p=sig*y2[i-1]+2.0; y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (opt) qn=un=0.0;
  else { qn=0.5; un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2])); }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2; k>=0; k--) y2[k]=y2[k]*y2[k+1]+u[k];

  /* Clean up */
  free(u);

  return 1;

}
