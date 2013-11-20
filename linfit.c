/*************************************************************************** 
LINFIT: This is the FIT algorithm taken from Numercial Recipes. Here's what
NR has to say about it

"Given a set of data points x[1..ndata],y[1..ndata] with individual
standard deviations sig[1..ndata], fit them to a straight line y = a + bx
by minimizing chisquared. Returned are a,b and their respective probable
uncertainties siga and sigb, chisquared chi2, and the goodness-of-fit
probability q (that the fit would have chisquared this large or larger). If
mwt=0 on input, then the standard deviations are assumed to be unavailable:
q is returned as 1.0 and the normalization of chi2 is to unit standard
deviation on all points."

I have modified the routine in the following ways:
0) Double precision and now an integer function
1) Added own error handling
2) Two options are now available:
    opt1=1 or 0 for when standard deviations on y do or do not exist respectively
    opt2=1 or 0 for when a is to vary or is to remain fixed (at input
         value) respectively
3) NOTE: When using opt2=0, a sum over x[i] is performed. Therefore,
   rounding errors could be important if x contains elements with widely
   varying values.
4) Added status array so that not all data points will be used. This
   can be passed as NULL is all points are to be used.
***************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "error.h"

int linfit(double *x, double *y, double *sig, int *sts, int ndata, int opt1, int opt2,
	   double *a, double *b, double *siga, double *sigb, double *chi2) {

  double wt=0.0,t=0.0,sxoss=0.0,sx=0.0,sy=0.0,st2=0.0,ss=0.0,sigdat=0.0;
  int    nval=0;
  int    i=0;

  /* Make sure number of points is reasonable */
  if (opt2 && ndata<2) {
    nferrormsg("linfit(): Cannot fit 2-parameter straight line\n\
\tto only %d data points",ndata); return 0;
  }
  if (!opt2 && ndata<1) {
    nferrormsg("linfit(): Cannot fit slope to only %d data points",ndata);
    return 0;
  }
  /* Find number of valid data points and make sure that's reasonable too */
  nval=ndata; if (sts!=NULL) {
    for (i=0; i<ndata; i++) { if (sts[i]!=1) nval--; }
  }
  if (opt2 && nval<2) {
    nferrormsg("linfit(): Cannot fit 2-parameter straight\n\
\tline to only %d valid data points",nval); return 0;
  }
  if (!opt2 && nval<1) {
    nferrormsg("linfit(): Cannot fit slope to only %d\n\
\tvalid data points",nval); return 0;
  }

  /* Determine slope and/or intercept */
  if (opt2) {
    *b=0.0;
    if (opt1) {
      ss=0.0; for (i=0; i<ndata; i++) {
	if (sts==NULL || sts[i]==1) {
	  if (!sig[i]) {
	    nferrormsg("linfit(): Sig array contains zero entry\n\tat index %d",i);
	    return 0;
	  }
	  wt=1.0/(sig[i]*sig[i]); ss+=wt; sx+=x[i]*wt; sy+=y[i]*wt;
	}
      }
    } else {
      for (i=0; i<ndata; i++) { if (sts==NULL || sts[i]==1) { sx+=x[i]; sy+=y[i]; } }
      ss=(double)nval;
    }
    sxoss=sx/ss;
    if (opt1) {
      for (i=0; i<ndata; i++) {
	if (sts==NULL || sts[i]==1) {
	  t=(x[i]-sxoss)/sig[i]; st2+=t*t; *b+=t*y[i]/sig[i];
	}
      }
    } else {
      for (i=0; i<ndata; i++) {
	if (sts==NULL || sts[i]==1) {
	  t=x[i]-sxoss; st2+=t*t; *b+=t*y[i];
	}
      }
    }
    *b/=st2; *a=(sy-sx*(*b))/ss; *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
    *sigb=sqrt(1.0/st2); *chi2=0.0;
    if (opt1) {
      for (i=0; i<ndata; i++) {
	if (sts==NULL || sts[i]==1)
	  *chi2+=((y[i]-(*a)-(*b)*x[i])/sig[i])*((y[i]-(*a)-(*b)*x[i])/sig[i]);
      }
    } else {
      for (i=0; i<ndata; i++) {
	if (sts==NULL || sts[i]==1) *chi2+=(y[i]-(*a)-(*b)*x[i])*(y[i]-(*a)-(*b)*x[i]);
      }
      *siga*=sigdat; *sigb*=sigdat;
    }
  } else {
    if (opt1) {
      for (i=0; i<ndata; i++) {
	if (sts==NULL || sts[i]==1) {
	  if (!sig[i]) {
	    nferrormsg("linfit(): Sig array contains zero entry\n\tat index %d",i);
	    return 0;
	  }
	  wt=1.0/(sig[i]*sig[i]); sx+=x[i]*x[i]*wt; sy+=x[i]*(y[i]-(*a))*wt;
	}
      }
    } else {
      for (i=0; i<ndata; i++) {
	if (sts==NULL || sts[i]==1) { sx=x[i]*x[i]; sy+=x[i]*(y[i]-(*a)); }
      }
    }
    *b=sy/sx; *sigb=sqrt(1.0/sx); *chi2=0.0;
    if (!opt1) {
      for (i=0; i<ndata; i++) {
	if (sts==NULL || sts[i]==1) *chi2+=(y[i]-(*a)-(*b)*x[i])*(y[i]-(*a)-(*b)*x[i]);
      }
      sigdat=sqrt((*chi2)/(nval-1)); *sigb*=sigdat;
    } else {
      for (i=0; i<ndata; i++) {
	if (sts==NULL || sts[i]==1)
	  *chi2+=((y[i]-(*a)-(*b)*x[i])/sig[i])*((y[i]-(*a)-(*b)*x[i])/sig[i]);
      }
    }
  }

  return 1;

}
