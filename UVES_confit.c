/****************************************************************************
* Find a continuum given data, err and clip-status arrays and a
* rejection significance for pixels above and below the continuum. The
* initial fit is made with the (pctl*100)th percentile removed and
* each subsequent iteration includes only data within the rejection
* limits of the current fit. Iterations continue until no more points
* are removed between iterations. If the err array is passed as NULL
* then an artificial err array is created which is simply unity
* everywhere. The sigma-clipping thresholds are ignored and the fit is
* only done once with all the data (i.e. no percentiles are removed).
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "fit.h"
#include "sort.h"
#include "memory.h"
#include "error.h"

#define FREE_ALL free(coeff); free(fit_x); free(fit_d); free(fit_e);

int UVES_confit(double *dat, double *err, int *sts, int ndat, int fit_typ,
		int fit_ord, double lrejsig, double urejsig, double pctl,
		int verb, double *fit) {

  double   pctl_min=0.0,chisq=0.0;
  double   *coeff=NULL,*fit_x=NULL,*fit_d=NULL,*fit_e=NULL;
  double   *work_a1=NULL,**work_m1=NULL,**work_m2=NULL;
  int      fit_n=0,fit_nr=0,iter=0,npctl=0;
  int      i=0;

  /* Make sure that the percentile range is reasonable */
  if (err!=NULL && (pctl<0.0 || pctl>=1.0)) {
    nferrormsg("UVES_confit(): Percentile range entered (=%lf)\n\
\tmust be >=0.0 and <1.0",pctl*100.0); return 0;
  }

  /* Allocate memory for fitting arrays */
  if ((coeff=darray(fit_ord))==NULL)
    errormsg("UVES_confit(): Could not allocate memory\n\
\tfor coeff array of size %d",fit_ord);
  if ((fit_x=darray(ndat))==NULL)
    errormsg("UVES_confit(): Could not allocate memory\n\
\tfor fit_x array of size %d",ndat);
  if ((fit_d=darray(ndat))==NULL)
    errormsg("UVES_confit(): Could not allocate memory\n\
\tfor fit_d array of size %d",ndat);
  if ((fit_e=darray(ndat))==NULL)
    errormsg("UVES_confit(): Could not allocate memory\n\
\tfor fit_e array of size %d",ndat);
  
  /* Set up for initial fit */
  for (fit_n=0,i=0; i<ndat; i++) if (sts[i]==1) fit_d[fit_n++]=dat[i];
  npctl=(err==NULL) ? fit_n : fit_n-(int)(pctl*(double)fit_n);
  if (npctl<=fit_ord) {
    if (verb) nferrormsg("UVES_confit(): Cannot make fit of order %d to\n\
\t%d VALID data points left after lowest %6.3lfth percentile removed",
			 fit_ord,npctl,pctl*100.0);
    FREE_ALL; return -2;
  }
  if (npctl!=fit_n) pctl_min=dselect(fit_n-npctl,fit_n,fit_d);
  else pctl_min=-INFIN;
  if (fit_typ==FITCHE) {
    for (fit_n=0,i=0; i<ndat; i++) {
      if (sts[i]==1 && dat[i]>=pctl_min) {
	fit_x[fit_n]=2.0*(double)i/(double)(ndat-1)-1.0; fit_d[fit_n]=dat[i];
	fit_e[fit_n++]=(err==NULL) ? 1.0 : err[i];
      }
    }
  } else {
    for (fit_n=0,i=0; i<ndat; i++) {
      if (sts[i]==1 && dat[i]>=pctl_min) {
	fit_x[fit_n]=(double)i; fit_d[fit_n]=dat[i];
	fit_e[fit_n++]=(err==NULL) ? 1.0 : err[i];
      }
    }
  }

  /* Now iterate a polynomial fit, rejecting points rejsig-sigma below the
     continuum */
  fit_nr=0; iter=0;
  while (fit_n!=fit_nr && fit_n>fit_ord && iter<MAXITERFIT) {
    fit_nr=fit_n;
    if (fit_typ==FITPOL) {
      if (!svdfit(fit_x,fit_d,fit_e,fit_n,coeff,fit_ord,&work_m1,&work_m2,&work_a1,
		  &chisq,svdfit_poly)) {
	nferrormsg("UVES_confit(): Error returned from svdfit()\n\
\twhile fitting regular polynomial"); return 0;
      }
      for (i=0; i<ndat; i++) fit[i]=poly_eval((double)i,coeff,fit_ord);
    } else if (fit_typ==FITCHE) {
      if (!svdfit(fit_x,fit_d,fit_e,fit_n,coeff,fit_ord,&work_m1,&work_m2,&work_a1,
		  &chisq,svdfit_chebyshev)) {
	nferrormsg("UVES_confit(): Error returned from svdfit()\n\
\twhile fitting Chebyshev polynomial"); return 0;
      }
      for (i=0; i<ndat; i++)
	fit[i]=chebyshev_eval(2.0*(double)i/(double)(ndat-1)-1.0,coeff,fit_ord);
    } else if (fit_typ==FITLEG) {
      if (!svdfit(fit_x,fit_d,fit_e,fit_n,coeff,fit_ord,&work_m1,&work_m2,&work_a1,
		  &chisq,svdfit_legendre)) {
	nferrormsg("UVES_confit(): Error returned from svdfit()\n\
\twhile fitting Legendre polynomial"); return 0;
      }
      for (i=0; i<ndat; i++) fit[i]=legendre_eval((double)i,coeff,fit_ord);
    }
    free(work_a1); free(*work_m1); free(work_m1); free(*work_m2); free(work_m2);
    if (fit_typ==FITCHE) {
      for (fit_n=0,i=0; i<ndat; i++) {
	if (sts[i]==1 && (err==NULL || (dat[i]>fit[i]-lrejsig*err[i] &&
					dat[i]<fit[i]+urejsig*err[i]))) {
	  fit_x[fit_n]=2.0*(double)i/(double)(ndat-1)-1.0; fit_d[fit_n]=dat[i];
	  fit_e[fit_n++]=(err==NULL) ? 1.0 : err[i];
	}
      }
    }
    else {
      for (fit_n=0,i=0; i<ndat; i++) {
	if (sts[i]==1 && (err==NULL || (dat[i]>fit[i]-lrejsig*err[i] &&
					dat[i]<fit[i]+urejsig*err[i]))) {
	  fit_x[fit_n]=(double)i; fit_d[fit_n]=dat[i];
	  fit_e[fit_n++]=(err=NULL) ? 1.0 : err[i];
	}
      }
    }
    iter++;
  }
  /* Check exit status from loop */
  if (fit_n<=fit_ord) {
    if (verb) nferrormsg("UVES_confit(): Only %d pixels left for continuum fit",fit_n);
    FREE_ALL; return -2;
  } else if (iter==MAXITERFIT) {
    if (verb) nferrormsg("UVES_confit(): Failed to converge on unique continuum\n\
\tin %d iterations. Change fit order or number of fitted points",MAXITERFIT);
    FREE_ALL; return -1;
  }

  /* Construct final continuum fit */
  if (fit_typ==FITPOL)
    for (i=0; i<ndat; i++) fit[i]=poly_eval((double)i,coeff,fit_ord);
  else if (fit_typ==FITCHE)
    for (i=0; i<ndat; i++)
      fit[i]=chebyshev_eval(2.0*(double)i/(double)(ndat-1)-1.0,coeff,fit_ord);
  else if (fit_typ==FITLEG)
    for (i=0; i<ndat; i++) fit[i]=legendre_eval((double)i,coeff,fit_ord);

  /* Clean up */
  FREE_ALL;
  
  return 1;

}
