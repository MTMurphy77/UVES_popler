/*************************************************************************** 
SVDVAR: This is the SVDVAR algorithm taken from Numercial Recipes. Here's
what NR has to say about it

To evaluate the covariance matrix cvm[1..ma][1..ma] of the fit for ma
parameters obtained by svdfit. Call this routine with matrices
v[1..ma][1..ma], w[1..ma] as returned from svdfit.

I have modified the routine in the following ways:

0) Double precision and now an integer function
1) Added own error handling

****************************************************************************/

#include <stdlib.h>
#include "memory.h"
#include "error.h"

int svdvar(double **v, int ma, double *w, double ***cvm) {

  int    i=0,j=0,k=0;
  double sum=0.0;
  double *wti=NULL;

  /* Allocate memory for relevant arrays */
  if ((wti=darray(ma))==NULL) {
    nferrormsg("svdvar(): Cannot allocate memory to wti\n\tarray of size %d",ma);
    return 0;
  }

  /* Allocate memory for covariance matrix */
  if ((*cvm=dmatrix(ma,ma))==NULL) {
    nferrormsg("svdvar(): Cannot allocate memory to matrix\n\t of size %dx%d",
	       ma,ma); return 0;
  }

  for (i=0; i<ma; i++) {
    wti[i]=0.0;
    if (w[i]) wti[i]=1.0/(w[i]*w[i]);
  }
  for (i=0; i<ma; i++) {
    for (j=0; j<=i; j++) {
      for (sum=0.0,k=0; k<ma; k++) sum+=v[i][k]*v[j][k]*wti[k];
      (*cvm)[j][i]=(*cvm)[i][j]=sum;
    }
  }

  /* Clean up */
  free(wti);

  return 1;

}
