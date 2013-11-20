/*************************************************************************** 
MRQCOF: This is the mrqcof algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Used by mrqmin to evaluate the linearized fitting matrix alpha, and vector
beta as in (15.5.8), and calculate chisq."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.
2) First argument must now be a matrix of values, x=x[ndata][nx]. The
   fifth argument, nx, specifies the number of x-points necessary to
   calculate the fit (and it's derivatives) in the function funcs. For
   example, for mrqfit_gauss(), nx=1 is all that's required. However,
   for fitting a Gaussian which has been sampled by pixels, as in
   mqrfit_erffn(), integration across the pixels is required and so
   two x-values are required per pixel (the start and ending
   wavelength of the pixel).
3) External function funcs() now must be an integer function and has an
   additional matrix argument specifying the externally defined limits for
   the parameters
****************************************************************************/

#include <stdlib.h>
#include "fit.h"
#include "memory.h"
#include "error.h"

int mrqcof(double **x, double *y, double *sig, int ndata, int nx, double *a,
	   double **lim, int *ia, int **inf, int ma, double **alpha, 
	   double *beta, double *chisq,
	   int (*funcs)(double *, double *, double **, int *, int **, double *,
			double *, int, int, int)) {

  double ymod=0.0,wt=0.0,sig2i=0.0,dy=0.0,*dyda=NULL;
  int i=0,j=0,k=0,l=0,m=0,mfit=0;

  if ((dyda=darray(ma))==NULL)
    errormsg("mrqcof(): Cannot allocate memory for dyda array\n\
\tof size %d",ma);
  
  for (j=0; j<ma; j++) if (ia[j]) mfit++;
  for (j=0; j<mfit; j++) {
    beta[j]=0.0; for (k=0; k<=j; k++) alpha[j][k]=0.0;
  }
  *chisq=0.0;
  for (i=0; i<ndata; i++) {
    if(!((*funcs)(x[i],a,lim,ia,inf,&ymod,dyda,nx,ma,i)))
      errormsg("mrqcof(): Unkown error returned from (*funcs)()");
    sig2i=1.0/(sig[i]*sig[i]);
    dy=y[i]-ymod;
    for (j=0,l=0; l<ma; l++) {
      if (ia[l]) {
	wt=dyda[l]*sig2i;
	for (k=0,m=0; m<=l; m++) {
	  if (ia[m]) alpha[j][k++]+=wt*dyda[m];	
	}
	beta[j++]+=dy*wt;
      }
    }
    *chisq+=dy*dy*sig2i;
  }
  for (j=1; j<mfit; j++) {
    for (k=0; k<j; k++) alpha[k][j]=alpha[j][k];
  }

  free(dyda);

  return 1;

}
