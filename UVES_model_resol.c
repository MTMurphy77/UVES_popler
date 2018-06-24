/****************************************************************************
* Use the ThAr line parameters and seeing information from the ThAr
* and object extractions to model the final resolution of UVES as a
* function of position across each echelle order. For each ThAr line
* used in the polynomial wavelength solution, the intrinsic
* instrumental line-spread function (ILSF) is estimated by subtracting
* (in quadrature) the equivalent Gaussian FWHM of the slit-width
* (projected on the CCD) from the measured ThAr line FWHM. A fit is
* then made to the LSF width versus position along the echelle
* order. From this fit, the LSF width for a given pixel in an echelle
* order is added in quadrature with the equivalent Gaussian FWHM of
* the seeing profile truncated by the slit. This gives a decent
* approximation to the expected spectral resolution as a fnuction of
* position along each echelle order.
* 
* Note that the resulting resolution arrays calculated here are still
* in pixel space. They should be converted to wavelength or velocity
* space once the dispersion for each pixel is determined from the
* wavelength solution.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "fit.h"
#include "stats.h"
#include "memory.h"
#include "const.h"
#include "error.h"

#define FREE_SVDWORK free(work_a1); free(*work_m1); free(work_m1); free(*work_m2); \
  free(work_m2);
#define FREE_FIT free(fit_d); free(fit_x); free(fit_e); free(sts);

int UVES_model_resol(spectrum *spec) {

  double  rejsig=3.0;
  double  pixscal=0.0,FWHM_SW=0.0,FWHM_SE=0.0,x=0.0,x2=0.0,x4=0.0,chisq=0.0,rms=0.0;
  double  *coeff=NULL,*fit_x=NULL,*fit_d=NULL,*fit_e=NULL;
  double  *wlsfsq=NULL;
  double  *work_a1=NULL,**work_m1=NULL,**work_m2=NULL;
  int     fit_ord=4,niter=3;
  int     fit_n=0;
  int     i=0,j=0,k=0;
  int     *sts=NULL;
  statset stat;

  /* Check that there's some useful ThAr lines to use */
  if (spec->ts.n<fit_ord || spec->ts.np<fit_ord) {
    nferrormsg("UVES_model_resol(): Number of ThAr lines in parameter\n\
\tlist (as read from from\n\t%s)\n\t(%d) or used in polynomial wavelength solution\n\
\t(%d) is too small. There may have been a problem reading the line\n\
\tparameters in the above file",spec->wlfile,spec->ts.n,spec->ts.np);
    return 0;
  }

  /* Loop over ThAr lines and estimate line-spread function widths */
  for (i=0,j=0; i<spec->ts.n; i++) {
    /* Only consider ThAr lines used in the polynomial solution */
    if (spec->ts.stp[i]) {
      /* Estimate width of intrinsic instrumental line-spread function by
	 subtracting (in quadrature) from the measured ThAr FWHM the
	 equivalent Gaussian FWHM of the slit-width projected to this
	 position on the CCD */
      pixscal=UVES_pixscal(spec->cwl,spec->ts.x[i],spec->or[0].np,spec->ts.binx);
      FWHM_SW=C_SQRT3*spec->wc_sw/pixscal/2.0;
      spec->ts.wlsf[i]=spec->ts.w[i]*spec->ts.w[i]-FWHM_SW*FWHM_SW;
    }
  }

  /* Allocate memory for fitting arrays */
  if ((coeff=darray(fit_ord))==NULL)
    errormsg("UVES_model_resol(): Could not allocate memory\n\
\tfor coeff array of size %d",fit_ord);
  if ((fit_d=darray(spec->ts.n))==NULL)
    errormsg("UVES_model_resol(): Could not allocate memory\n\
\tfor fit_d array of size %d",spec->ts.n);
  if ((fit_x=darray(spec->ts.np))==NULL)
    errormsg("UVES_model_resol(): Could not allocate memory\n\
\tfor fit_x array of size %d",spec->ts.np);
  if ((fit_e=darray(spec->ts.np))==NULL)
    errormsg("UVES_model_resol(): Could not allocate memory\n\
\tfor fit_e array of size %d",spec->ts.np);
  if ((sts=iarray(spec->ts.n))==NULL)
    errormsg("UVES_model_resol(): Cannot allocate memory for sts\n\
\tarray of length %d",spec->ts.n);

  /* Create initial constant error array by taking the RMS of all LSF measurements */
  if (!stats(spec->ts.wlsf,NULL,NULL,NULL,spec->ts.stp,spec->ts.n,0,&stat)) {
    nferrormsg("UVES_model_resol(): Error returned from stats()\n\
\twhen attempting to find RMS of LSF measurements from file\n\t%s",spec->wlfile);
    return 0;
  }
  /* Fill error array and status array for use in sigma-clipping algorithm */
  for (i=0; i<spec->ts.np; i++) fit_e[i]=stat.rms;
  for (i=0; i<spec->ts.n; i++) sts[i]=(spec->ts.stp[i]) ? 1 : 0;
  /* Iterate the fit */
  fit_n=spec->ts.np;
  for (j=0; j<niter; j++) {
    /* Fill fitting arrays according to current status array */
    for (i=0,k=0; i<spec->ts.n; i++) {
      if (spec->ts.stp[i] && sts[i]) {
	fit_x[k]=2.0*(spec->ts.x[i]-1.0)/(double)(spec->or[0].np-1)-1.0;
	fit_d[k++]=spec->ts.wlsf[i];
      }
    }
    /* Do the fit */
    if (!svdfit(fit_x,fit_d,fit_e,fit_n,coeff,fit_ord,&work_m1,&work_m2,&work_a1,
		&chisq,svdfit_chebyshev)) {
      nferrormsg("UVES_model_resol(): Error returned from svdfit()\n\
\twhile fitting Chebyshev polynomial"); return 0;
    }
    for (i=0; i<spec->ts.n; i++)
      fit_d[i]=chebyshev_eval(2.0*(spec->ts.x[i]-1.0)/(double)(spec->or[0].np-1)-1.0,
			      coeff,fit_ord);
    FREE_SVDWORK;
    if (j==niter-1) break;
    /* Calculate RMS of valid points around current fit */
    for (i=0,rms=0.0; i<spec->ts.n; i++) if (spec->ts.stp[i] && sts[i])
      rms+=(spec->ts.wlsf[i]-fit_d[i])*(spec->ts.wlsf[i]-fit_d[i]);
    rms/=(double)fit_n; rms=sqrt(rms);
    /* Reject deviant points, calculate new number of valid points and
       reset error array */
    for (i=0,fit_n=0; i<spec->ts.n; i++) {
      if (spec->ts.stp[i] && fabs(spec->ts.wlsf[i]-fit_d[i])<rejsig*rms) {
	sts[i]=1; fit_e[fit_n++]=rms;
      } else sts[i]=0;
    }
    /* Check that number of valid points is reasonable */
    if (fit_n<fit_ord) {
      nferrormsg("UVES_model_resol(): Number of valid data points\n\
\t(=%d out of a total of %d) is too low for a fit of order %d\n\
\tafter iteration %d",fit_n,spec->ts.np,j+1);
      FREE_FIT; return 0;
    }
  }
  /* Clean up */
  FREE_FIT;

  /* Using the fit coefficients, construct a LSF width array of the
     same dimensions as the echelle orders */
  /* Allocate memory */
  if ((wlsfsq=darray(spec->or[0].np))==NULL)
    errormsg("UVES_model_resol(): Could not allocate memory\n\
\tfor wlsfsq array of size %d",spec->or[0].np);
  /* Construct LSF width array */
  for (i=0; i<spec->or[0].np; i++)
    wlsfsq[i]=chebyshev_eval(2.0*(double)i/(double)(spec->or[0].np-1)-1.0,coeff,
			     fit_ord);

  /* Construct resolution array for each echelle order. The LSF width
     is added in quadrature to the equivalent Gaussian width of the
     seeing profile truncated by the slit as projected to the relevent
     position along the order */
  for (i=0; i<spec->nor; i++) {
    if (spec->or[i].nuse>=MINUSE) {
      x=spec->or[i].seeing/spec->sw; x2=x*x; x4=x2*x2;
      x=1.0/pow((1.0/x4+16.0/9.0),0.25)+1.05/x4/x2/(exp(4.0/x)-1.0)-
	3.3/x4/(exp(4.4/x)-1.0);
      for (j=0; j<spec->or[i].np; j++) {
	pixscal=UVES_pixscal(spec->cwl,(double)(j+1),spec->or[i].np,spec->binx);
	FWHM_SE=x*spec->sw/pixscal;
	spec->or[i].res[j]=sqrt(wlsfsq[j]+FWHM_SE*FWHM_SE);
	spec->or[i].res[j]=(MAX(spec->or[i].res[j],0.0));
      }
    }
  }
  /* Clean up */
  free(wlsfsq); free(coeff);
  
  return 1;

}
