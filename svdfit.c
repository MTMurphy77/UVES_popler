/*************************************************************************** 
SVDFIT: This is the SVDFIT algorithm taken from Numercial Recipes. Here's
what NR has to say about it

Given a set of data points x[1..ndata],y[1..ndata] with individual standard
deviations sig[1..ndata], use chisquared minimization to determine the
coefficients a[1..ma] of the fitting function y = sum{i} a_i ×
afunci(x). Here we solve the fitting equations using singular value
decomposition of the ndata x ma matrix, as in §2.6. Arrays
u[1..ndata][1..ma], v[1..ma][1..ma], and w[1..ma] provide workspace on
input; on output they define the singular value decomposition, and can be
used to obtain the covariance matrix. The program returns values for the ma
fit parameters a, and chisq. The user supplies a routine funcs(x,afunc,ma)
that returns the ma basis functions evaluated at x = x in the array
afunc[1..ma].

I have modified the routine in the following ways:

0) Double precision and now an integer function
1) Made the function "funcs" an integer function
2) Memory allocated to u, v and w instead of assuming they are passed
   already filled
3) Added own error handling
****************************************************************************/

#include <stdlib.h>
#include "fit.h"
#include "memory.h"
#include "error.h"

int svdfit(double *x, double *y, double *sig, int ndata, double *a, int ma,
	double ***u, double ***v, double **w, double *chisq,
	int (*funcs)(double, double *, int)) {

  int     i=0,j=0;
  double  wmax=0.0,tmp=0.0,thresh=0.0,sum=0.0;
  double  *b=NULL,*afunc=NULL;

  /* Allocate memory for relevant arrays/matrices */
  if ((b=darray(ndata))==NULL) {
    nferrormsg("svdfit(): Cannot allocate memory to b\n\tarray of size %d",
	       ndata); return 0;
  }
  if ((afunc=darray(ma))==NULL) {
    nferrormsg("svdfit(): Cannot allocate memory to afunc\n\tarray of size %d",
	       ma); return 0;
  }
  if ((*w=darray(ma))==NULL) {
    nferrormsg("svdfit(): Cannot allocate memory to w\n\tarray of size %d",
	       ma); return 0;
  }
  if ((*u=dmatrix(ndata,ma))==NULL) {
    nferrormsg("svdfit(): Cannot allocate memory to matrix\n\t of size %dx%d",
	       ndata,ma); return 0;
  }
  if ((*v=dmatrix(ma,ma))==NULL) {
    nferrormsg("svdfit(): Cannot allocate memory to matrix\n\t of size %dx%d",
	       ma,ma); return 0;
  }

  /* Begin SVD fitting */
  for (i=0; i<ndata; i++) {
    if (!(*funcs)(x[i],afunc,ma)) {
      nferrormsg("svdfit(): Error returned from fitting function");
      return 0;
    }
    tmp=1.0/sig[i];
    for (j=0; j<ma; j++) (*u)[i][j]=afunc[j]*tmp;
    b[i]=y[i]*tmp;
  }
  if (!svdcmp(*u,ndata,ma,*w,*v)) {
    nferrormsg("svdfit(): Error returned from svdcmp()"); return 0;
  }
  wmax=0.0; for (j=0; j<ma; j++) wmax=MAX(wmax,(*w)[j]); thresh=SVDFIT_TOL*wmax;
  for (j=0; j<ma; j++) {
    if ((*w)[j]<thresh) {
      (*w)[j]=0.0;
      warnmsg("svdfit(): Setting coefficient %d's singular value to zero",j);
    }
  }
  if (!svbksb(*u,*w,*v,ndata,ma,b,a)) {
    nferrormsg("svdfit(): Error returned from svbksb()"); return 0;
  }
  *chisq=0.0;
  for (i=0; i<ndata; i++) {
    if (!(*funcs)(x[i],afunc,ma)) {
      nferrormsg("svdfit(): Error returned from fitting function");
      return 0;
    }
    for (sum=0.0,j=0; j<ma; j++) sum+=a[j]*afunc[j];
    *chisq+=(tmp=(y[i]-sum)/sig[i],tmp*tmp);
  }

  /* Clean up */
  free(b); free(afunc);

  return 1;

}
