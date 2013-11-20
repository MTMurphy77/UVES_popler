/*************************************************************************** 
ERFFN: This is the erff algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Returns the error function erf(x)."

I have changed the raw routine to include the following features:
0) Double precision is used.

****************************************************************************/

#include "gamm.h"

double erffn(double x) {

  return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);

}
