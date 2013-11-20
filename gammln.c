/****************************************************************************
GAMMLN: This is the logarithmic gamma function algorithm taken from
Numerical Recipes. Here's what NR has to say;

"Returns the value ln|gamma(x)| for x > 0.".

I have changed the raw routine to include the following:
0) Double precision is used

****************************************************************************/

#include <math.h>
#include "gamm.h"

double gammln(double x) {

  double        xx,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=xx=x;
  tmp=xx+5.5;
  tmp-=(xx+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0; j<=5; j++) ser+=cof[j]/++y;

  return -tmp+log(2.5066282746310005*ser/xx);
}
