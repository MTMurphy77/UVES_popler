/*************************************************************************** 
GAUSSJ: This is the gaussj algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Linear equation solution by Gauss-Jordan elimination, equation (2.1.1)
above. a[1..n][1..n] is the input matrix. b[1..n][1..m] is input containing
the m right-hand side vectors. On output, a is replaced by its matrix
inverse, and b is replaced by the corresponding set of solution vectors."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.

BJZ (26/5/2006): Have prevented this routine from forming singular
                 matrices to stop problem with mrqfit where solution
                 is degenerate. When this happens the routine returns
                 GJSNGLR rather than 0 (normal error) or 1 (normal
                 exit status)
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "fit.h"
#include "memory.h"
#include "error.h"

int gaussj(double **a, int n, double **b, int m) {

  double big=0.0,dum=0.0,pivinv=0.0,tempr=0.0;
  int    ret=1; /* Default return value */
  int    *indxc=NULL,*indxr=NULL,*ipiv=NULL;
  int    i=0,icol=0,irow=0,j=0,k=0,l=0,ll=0;

  if ((indxc=iarray(n))==NULL)
    errormsg("gaussj(): Cannot allocate memory for indxc array\n\
\tof size %d",n);
  if ((indxr=iarray(n))==NULL)
    errormsg("gaussj(): Cannot allocate memory for indxr array\n\
\tof size %d",n);
  if ((ipiv=iarray(n))==NULL)
    errormsg("gaussj(): Cannot allocate memory for ipiv array\n\
\tof size %d",n);

  for (i=0; i<n; i++) {
    big=0.0;
    for (j=0; j<n; j++) {
      if (ipiv[j]!=1) {
	for (k=0; k<n; k++) {
	  if (ipiv[k]==0) {
	    if ((dum=fabs(a[j][k]))>=big) { big=dum; irow=j; icol=k; }
	  }
	  else if (ipiv[k]>1) errormsg("gaussj(): Singular matrix (1) formed");
	}
      }
    }
    ++(ipiv[icol]);
    if (irow!=icol) {
      for (l=0; l<n; l++) { SWAP((a[irow][l]),(a[icol][l])); }
      for (l=0; l<m; l++) { SWAP((b[irow][l]),(b[icol][l])); }
    }
    indxr[i]=irow; indxc[i]=icol;
    if (a[icol][icol]==0.0) {
      a[icol][icol]=FPMIN; ret=GJSNGLR;
      nferrormsg("gaussj(): Singular matrix (2) formed - setting [%d][%d]\n\
\tentry to FPMIN to prevent division by zero.",icol,icol);
    }
    pivinv=1.0/a[icol][icol]; a[icol][icol]=1.0;
    for (l=0; l<n; l++) a[icol][l]*=pivinv;
    for (l=0; l<m; l++) b[icol][l]*=pivinv;
    for (ll=0; ll<n; ll++) {
      if (ll!=icol) {
	dum=a[ll][icol]; a[ll][icol]=0.0;
	for (l=0; l<n; l++) a[ll][l]-=a[icol][l]*dum;
	for (l=0; l<m; l++) b[ll][l]-=b[icol][l]*dum;
      }
    }
  }
  for (l=n-1; l>=0; l--) {
    if (indxr[l]!=indxc[l]) {
      for (k=0; k<n; k++) { SWAP((a[k][indxr[l]]),(a[k][indxc[l]])); }
    }
  }

  free(ipiv); free(indxr); free(indxc);

  return ret;

}
