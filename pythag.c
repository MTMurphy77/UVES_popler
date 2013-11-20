/*************************************************************************** 
PYTHAG: This is the PYTHAG algorithm taken from Numercial Recipes.

I have modified the routine in the following ways:

0) Double precision

****************************************************************************/

#include <math.h>
#include "fit.h"

double pythag(double a, double b) {

  double absa=0.0,absb=0.0;

  absa=fabs(a); absb=fabs(b);
  if (absa>absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));

}
