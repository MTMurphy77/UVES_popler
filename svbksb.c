/*************************************************************************** 
SVBKSB: This is the SVBKSB algorithm taken from Numercial Recipes. Here's
what NR has to say about it

Solves A.X = B for a vector X, where A is specified by the arrays
u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by svdcmp(). m and n are
the dimensions of a, and will be equal for square matrices. b[1..m] is the
input right-hand side. x[1..n] is the output solution vector. No input
quantities are destroyed, so the routine may be called sequentially with
different bs.

I have modified the routine in the following ways:

0) Double precision and now an integer function
1) Added own error handling

****************************************************************************/

#include <stdlib.h>
#include "memory.h"
#include "error.h"

int svbksb(double **u, double *w, double **v, int m, int n, double *b,
	   double *x) {

  double s=0.0;
  double *tmp=NULL;
  int    i=0,j=0,jj=0;
  
  /* Allocate memory for relevant arrays */
  if ((tmp=darray(n))==NULL) {
    nferrormsg("svbksb(): Cannot allocate memory to rv1\n\tarray of size %d",
	       n); return 0;
  }

  for (j=0; j<n; j++) {
    s=0.0;
    if (w[j]) {
      for (i=0; i<m; i++) {
	s+=u[i][j]*b[i];
      }
      s/=w[j];
    }
    tmp[j]=s;
  }
  for (j=0; j<n; j++) {
    s=0.0;
    for (jj=0; jj<n; jj++) s+=v[j][jj]*tmp[jj];
    x[j]=s;
  }

  /* Clean up */
  free(tmp);
  
  return 1;

}
