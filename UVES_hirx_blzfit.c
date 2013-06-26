/****************************************************************************
Default is to fit a Chebyshev polynomial, of order fit_ord, to the
blaze data for HIREDUX files. A new status array is formed which is 1
only when both the object flux status array is 1 and when the blaze
flux is positive. Before the first fit is done, pixels with fluxes
below 1/initrej or above initrej times the median blaze value are
rejected. The fit is then done with iterative rejection: at each
iteration pixels with fluxes lying more than rejsig times the RMS
residual (between the blaze data and current fit) away from the fit
are not considered in the next iteration. niter iterations are
performed. Note that the actual order of the fit is determined by the
proportion of valid blaze data. Lower order fits are conducted when a
smaller proportion of the data is valid.

The alternative is to calculate a running median in a window of width
"order" pixels. The same iterative rejection algorithm as above is
used. When using this method a fairly low niter value should be used
(probably 2 or 3 is just fine). In this method the initial rejection
of pixels (below 1/initrej or above rejsig times the median) is not
performed.

NOTE: When using SVDFIT, an artificial constant error array is
used. This is set to unity. If the real errors are very different to
unity then there may be instabilities in the fits. To alleviate the
possible effect of the above, the blaze data are first divided by the
median blaze value. The fit returned by this routine is, by default,
normalized by this median value. If the original scale is required,
opt1 should be set to 1.

opt1 : Set to 0 to return normalized blaze fit; 1 to use original
scale (see note above).

opt2 : Set to 0 for a Chebyshev polynomial fit to the blaze data; 1 to
conduct a running median filter instead.

****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "fit.h"
#include "stats.h"
#include "memory.h"
#include "error.h"

#define FREE_ALL free(blztmp); free(coeff); free(fit_x); free(fit_d); free(fit_e); \
  free(sts);

int UVES_hirx_blzfit(echorder *ord, double *blz, double initrej, double rejsig,
		     int order, int niter, int opt1, int opt2, double *fit) {

  double   chisq=0.0,rms=0.0,resid=0.0;
  double   *blztmp=NULL,*coeff=NULL,*fit_x=NULL,*fit_d=NULL,*fit_e=NULL;
  double   *work_a1=NULL,**work_m1=NULL,**work_m2=NULL;
  int      fit_n=0,fit_ord=0;
  int      i=0,j=0,k=0;
  int      *sts;
  statset  stat;

  /* Allocate memory for status array */
  if ((sts=iarray(ord->np))==NULL)
    errormsg("UVES_hirx_blzfit(): Could not allocate memory\n\
\tfor sts array of size %d",ord->np);
  
  /* Determine number of valid pixels and fill status array */
  for (i=0,fit_n=0; i<ord->np; i++) {
    if (ord->st[i]==1 && blz[i]>0.0) { sts[i]=1; fit_n++; }
    else sts[i]=0;
  }

  /* Make sure number of valid pixels is usable */
  if (fit_n<2) {
    nferrormsg("UVES_hirx_blzfit(): Number of valid pixels (=%d\n\
\tout of a total of %d) is too low for a reasonable fit",fit_n,ord->np);
    free(sts); return 0;
  }

  /* Determine median blaze value */
  if (!median(blz,sts,ord->np,&stat,0)) {
    nferrormsg("UVES_hirx_blzfit(): Error returned from median()");
    free(sts); return 0;
  }
  if (stat.med<=0.0) stat.med=1.0;

  /* Determine actual order of fit to use. Algorithm is just to use
     the sqrt of the ratio of the number of valid pixels and the total
     number of pixels. This is just arbitrary but seems to work
     well. */
  fit_ord=(MAX(2,(NINT(order*(sqrt((double)fit_n/(double)ord->np))))));

  /* Allocate memory for fitting arrays */
  if ((blztmp=darray(ord->np))==NULL)
    errormsg("UVES_hirx_blzfit(): Could not allocate memory\n\
\tfor blztmp array of size %d",ord->np);
  if ((coeff=darray(fit_ord))==NULL)
    errormsg("UVES_hirx_blzfit(): Could not allocate memory\n\
\tfor coeff array of size %d",fit_ord);
  if ((fit_x=darray(fit_n))==NULL)
    errormsg("UVES_hirx_blzfit(): Could not allocate memory\n\
\tfor fit_x array of size %d",fit_n);
  if ((fit_d=darray(fit_n))==NULL)
    errormsg("UVES_hirx_blzfit(): Could not allocate memory\n\
\tfor fit_d array of size %d",fit_n);
  if ((fit_e=darray(fit_n))==NULL)
    errormsg("UVES_hirx_blzfit(): Could not allocate memory\n\
\tfor fit_e array of size %d",fit_n);

  /* Normalize input blz array */
  for (i=0; i<ord->np; i++) blztmp[i]=blz[i]/stat.med;

  /* Set the artificial error array to unity */
  for (i=0; i<fit_n; i++) fit_e[i]=1.0;

  /* Identify the region over which a fit should be conducted and
     reject points which lie far away from median value */
  if (niter>0 && !opt2) {
    for (i=0,fit_n=0; i<ord->np; i++) {
      if (ord->st[i]==1 && blztmp[i]>0.0 && blztmp[i]>1.0/initrej &&
	  blztmp[i]<initrej) {
	sts[i]=1; fit_n++;
      } else sts[i]=0;
    }
  }

  /* Iterate the fit */
  for (j=0; j<niter; j++) {
    /* Do the fit or the running median */
    if (!opt2) {
      /* Fill fitting arrays according to current status array */
      for (i=0,k=0; i<ord->np; i++) {
	if (sts[i]) {
	  fit_x[k]=2.0*(double)i/(double)(ord->np-1)-1.0; fit_d[k++]=blztmp[i];
	}
      }
      /* Do the fit */
      if (!svdfit(fit_x,fit_d,fit_e,fit_n,coeff,fit_ord,&work_m1,&work_m2,&work_a1,
		  &chisq,svdfit_chebyshev)) {
	nferrormsg("UVES_hirx_blzfit(): Error returned from svdfit()\n\
\twhile fitting Chebyshev polynomial"); return 0;
      }
      for (i=0; i<ord->np; i++)
	fit[i]=chebyshev_eval(2.0*(double)i/(double)(ord->np-1)-1.0,coeff,fit_ord);
      free(work_a1); free(*work_m1); free(work_m1); free(*work_m2); free(work_m2);
    } else {
      /* Do the running median */
      if (!medianrun(blztmp,fit,sts,ord->np,order)) {
	  nferrormsg("UVES_hirx_blzfit(): Error returned from medianrun()\n\
\twhile calculating running median"); return 0;
      }
    }
    if (j==niter-1) break;
    /* Calculate RMS of valid points around current fit */
    for (i=0,rms=0.0; i<ord->np; i++) {
      if (sts[i]) { resid=blztmp[i]-fit[i]; rms+=resid*resid; }
    }
    rms/=(double)fit_n; rms=sqrt(rms);
    /* Reject deviant points and calculate new number of valid points */
    for (i=0,fit_n=0; i<ord->np; i++) {
      if (ord->st[i]==1 && blztmp[i]>0.0 && fabs(blztmp[i]-fit[i])<rejsig*rms) {
	sts[i]=1; fit_n++;
      } else sts[i]=0;
    }
    /* Check that number of valid points is reasonable */
    if (fit_n<fit_ord) {
      nferrormsg("UVES_hirx_blzfit(): Number of valid pixels (=%d\n\
\tout of a total of %d) is too low for a fit of order (or\n\
\trunning median of width) order %d after iteration %d",fit_n,ord->np,j+1);
      FREE_ALL; return 0;
    }
  }

  /* If required, renormalize the fit back to the original scale */
  if (opt1==1) for (i=0; i<ord->np; i++) fit[i]*=stat.med;

  /* Clean up */
  FREE_ALL;
  
  return 1;

}
