/*************************************************************************** 
SPLINT: This is the splint algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with
the x-axis in order), and given the array y2a[1..n], which is the output
from spline above, and given a value of x, this routine returns a
cubic-spline interpolated value y."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns a double, the interpolated value y.
2) Arrays are 0-indexed
3) Error checking done with own error codes
****************************************************************************/

#include "error.h"

double splint(double *xa, double *ya, double *y2a, int n, double x) {

  double h=0.0,b=0.0,a=0.0,y=0.0;
  int    klo=0,khi=0,k=0;

  /* Check to make sure array length is reasonable */
  if (n<2) errormsg("splint(): Number of data points (=%d) must\n\
\tbe greater than 1",n);

  klo=0; khi=n-1;
  while (khi-klo>1) {
    k=(khi+klo) >> 1;
    if (xa[k]>x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h<=0.0) errormsg("splint(): Input x array not in increasing order");
  a=(xa[khi]-x)/h; b=(x-xa[klo])/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;

  return y;

}
