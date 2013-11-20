/*************************************************************************** 
MRQMIN: This is the mrqmin algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Levenberg-Marquardt method, attempting to reduce the value chisq. of
a fit between a set of data points x[1..ndata], y[1..ndata] with
individual standard deviations sig[1..ndata], and a nonlinear function
dependent on ma coefficients a[1..ma]. The input array ia[1..ma]
indicates by nonzero entries those components of a that should be
fitted for, and by zero entries those components that should be held
fixed at their input values. The program returns current best-fit
values for the parameters a[1..ma], and chisq = chisq. The arrays
covar[1..ma][1..ma], alpha[1..ma][1..ma] are used as working space
during most iterations. Supply a routine funcs(x,a,yfit,dyda,ma) that
evaluates the fitting function yfit, and its derivatives dyda[1..ma]
with respect to the fitting parameters a, at x. On the first call
provide an initial guess for the parameters a, and set alamda<0 for
initialization (which then sets alamda=.001). If a step succeeds chisq
becomes smaller and alamda decreases by a factor of 10. If a step
fails alamda grows by a factor of 10. You must call this routine
repeatedly until convergence is achieved. Then, make one final call
with alamda=0, so that covar[1..ma][1..ma] returns the covariance
matrix, and alpha the curvature matrix. (Parameters held fixed will
return zero covariances.)"

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
3) "Spheres of influence" may now be specified for each
   parameter. This information is passed via a matrix, inf[1..ma][2],
   which contains the first and last pixels for which each parameter
   is relevant. E.g. in this manner one could specify a broken power
   law to fit the data. This matrix is passed directly to the funcs()
   routine where the information is interpretted appropriately for the
   given function.
4) External function funcs() now must be an integer function and has
   an additional matrices arguments specifying the externally defined
   limits for the parameters and their spheres of influence on the
   data.
5) gaussj() routine now returns GJSNGLR if there is a singular matrix
   formed (i.e. we have a degenerate solution). This same value is
   returned out of mrqmin if it occurs.
****************************************************************************/

#include <stdlib.h>
#include "fit.h"
#include "memory.h"
#include "error.h"

int mrqmin(double **x, double *y, double *sig, int ndata, int nx, double *a,
	   double **lim, int *ia, int **inf, int ma, double **covar,
	   double **alpha, double *chisq, double *alamda,
	   int (*funcs)(double *, double *, double **, int *, int **, double *,
			double *, int, int, int)) {

  static double ochisq,*atry,*beta,*da,**oneda;
  static int mfit;
  int ret=1; /* Default return value */
  int j=0,k=0,l=0;

  if (*alamda<0.0) {
    if ((atry=darray(ma))==NULL)
      errormsg("mrqmin(): Cannot allocate memory for atry array\n\
\tof size %d",ma);
    if ((beta=darray(ma))==NULL)
      errormsg("mrqmin(): Cannot allocate memory for beta array\n\
\tof size %d",ma);
    if ((da=darray(ma))==NULL)
      errormsg("mrqmin(): Cannot allocate memory for da array\n\
\tof size %d",ma);

    for (mfit=0,j=0; j<ma; j++) if (ia[j]) mfit++;
    if ((oneda=dmatrix(mfit,1))==NULL)
      errormsg("mrqmin(): Cannot allocate memeory for oneda matrix\n\
\tof size %dx%d",mfit,1);
    *alamda=0.001;
    if(!mrqcof(x,y,sig,ndata,nx,a,lim,ia,inf,ma,alpha,beta,chisq,funcs))
      errormsg("mrqmin(): Unknown error returned from mrqcof");
    ochisq=(*chisq);
    for (j=0; j<ma; j++) atry[j]=a[j];
  }

  for (j=0; j<mfit; j++) {
    for (k=0; k<mfit; k++) covar[j][k]=alpha[j][k];
    covar[j][j]=alpha[j][j]*(1.0+(*alamda)); oneda[j][0]=beta[j];
  }

  if (gaussj(covar,mfit,oneda,1)==GJSNGLR) ret=GJSNGLR;

  for (j=0; j<mfit; j++) da[j]=oneda[j][0];
  if (*alamda==0.0) {
    covsrt(covar,ma,ia,mfit); covsrt(alpha,ma,ia,mfit);
    free(*oneda); free(oneda); free(da); free(beta); free(atry);
    return 1;
  }

  for (j=0,l=0; l<ma; l++) if (ia[l]) atry[l]=a[l]+da[j++];
  mrqcof(x,y,sig,ndata,nx,atry,lim,ia,inf,ma,covar,da,chisq,funcs);

  if (*chisq < ochisq) {
    (*alamda)*=0.1; ochisq=(*chisq);
    for (j=0; j<mfit; j++) {
      beta[j]=da[j]; for (k=0; k<mfit; k++) alpha[j][k]=covar[j][k];
    }
    for (l=0; l<ma; l++) a[l]=atry[l];
  }
  else { (*alamda)*=10.0; *chisq=ochisq; }

  return ret;

}
