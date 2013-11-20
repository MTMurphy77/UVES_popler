/*************************************************************************** 
COVSRT: This is the covsrt algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Expand in storage the covariance matrix covar, so as to take into account
parameters that are being held fixed. (For the latter, return zero
covariances.)"

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.
****************************************************************************/

#include "fit.h"

int covsrt(double **covar, int ma, int *ia, int mfit) {

  double tempr=0.0;
  int i=0,j=0,k=0;
  
  for (i=mfit; i<ma; i++) {
    for (j=0; j<=i; j++) covar[i][j]=covar[j][i]=0.0;
  }
  for (k=mfit-1,j=ma-1; j>=0; j--) {
    if (ia[j]) {
      for (i=0; i<ma; i++) { SWAP((covar[i][k]),(covar[i][j])); }
      for (i=0; i<ma; i++) { SWAP((covar[k][i]),(covar[j][i])); }
      k--;
    }
  }

  return 1;

}
