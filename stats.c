/****************************************************************************
STATS: Function to calculate various statistics for given input arrays.

DESCRIPTION: The user must input a data array, dat, of length ndat. All
other double precision arrays can be left as NULL if the user doesn't
wish to calculate any weighted quantities. If the status array, sts,
is not NULL then it should be 1 when the pixel is valid and any other
quantity when invalid. If it is NULL then all pixels are assumed to be
valid. If the user enters a non-NULL sigma array, sig, but keeps the
expected fluctuation array, efl, and weights array, wgt, as NULL then
weighted quantites are calculated using 1/sigma^2 weighting and
chisq. is calculated assuming that sigma also estimates the expected
fluctuations. If the efl or wgt arrays are non-NULL then they are used
in the appropriate ways. The weighting used is 1/wgt^2 (i.e. the wgt
array is really the inverse square-root of the weights really used).

USAGE: opt = 0 when only the mean, rms, positive and negative rms and
               error in the mean are to be calculated (i.e. unweighted
               quantities).
       opt = 1 when all weighted quantities are to be calculated as well.
       opt = 2 when the median (and semi-interquartile range) is to be
               calculated in addition to weighted quantities.
       opt = 3 when the median sigma is to be calculated in addition
               to all the above.

****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "stats.h"
#include "error.h"

int stats(double *dat, double *sig, double *efl, double *wgt, int *sts, int ndat,
	  int opt, statset *stat) {

  double       sum=0.0,sumdatonsigsq=0.0,sumeflonsigsq=0.0,sumsigonsigsq=0.0;
  double       sum1onsigsq=0.0;
  double       dummy=0.0,dumsq=0.0;
  int          snull=0,nval=0,npos=0,nneg=0;
  register int i=0;

  /** Do for all values of opt **/
  /* Check that data array is usable */
  if (dat==NULL) { nferrormsg("stats(): Data array passed as NULL"); return 0; }

  /* Is STS to be used? */
  if (sts==NULL) snull=1;

  /* Initialise unweighted quantities */
  stat->rms=stat->prms=stat->nrms=stat->mean=stat->emean=0.0;

  /* Calculate mean of data */
  for (i=0; i<ndat; i++) if (snull || sts[i]==1) { sum+=dat[i]; nval++; }
  if (nval) stat->mean=sum/((double)nval);

  /* Calculate RMS, positive RMS and negative RMS of data */
  if (nval) {
    for (i=0; i<ndat; i++) {
      if (snull || sts[i]==1) {
	dummy=dat[i]-stat->mean; dumsq=dummy*dummy; stat->rms+=dumsq;
	if (dummy<0.0) { stat->nrms+=dumsq; nneg++; }
	else { stat->prms+=dumsq; npos++; }
      }
    }
    stat->rms=(nval>1) ? sqrt(stat->rms/((double)(nval-1))) : 0.0;
    if (!npos) stat->prms=0.0;
    else if (npos==1) stat->prms=sqrt(stat->prms);
    else stat->prms=sqrt(stat->prms/((double)(npos-1)));
    if (!nneg) stat->nrms=0.0;
    else if (nneg==1) stat->nrms=sqrt(stat->nrms);
    else stat->nrms=sqrt(stat->nrms/((double)(nneg-1)));
    /* Calculate error in mean of data */
    stat->emean=stat->rms/sqrt((double)nval);
  }

  /* Return if we don't want to calculate weighted quantities */
  if (!opt) return 1;

  /** Do if we want to calculate weighted quantities **/
  stat->meansig=stat->rmssig=stat->wmean=stat->ewmean=stat->eflwmean=0.0;
  stat->chisq=stat->rchisq=0.0;

  /* Return here if nval is zero */
  if (!nval && opt==1) return 1;
  else if (!nval && opt==2) { stat->med=stat->siqr=0.0; return 1; }
  else if (!nval && opt==3) { stat->med=stat->siqr=stat->emed=0.0; return 1; }

  /* Check to make sure we have an non-NULL sig array */
  if (sig==NULL) {
    nferrormsg("stats(): Sigma array passed as NULL.\n\
Must be non-NULL when using opt = %d\t",opt); return 0;
  }

  /* Calculate mean 1-sigma error, weighted mean, error in weighted
     mean and expected fluctuation in weighted mean */
  /** Diverge here for use of weights **/
  if (wgt==NULL) {
    /** Use sigma arrays as weights **/
    for (i=0,sum=0.0; i<ndat; i++) {
      if (snull || sts[i]==1) {
	sum+=sig[i]; dummy=sig[i]*sig[i]; sumdatonsigsq+=dat[i]/dummy;
	if (efl!=NULL) sumeflonsigsq+=efl[i]*efl[i]/dummy/dummy;
	sum1onsigsq+=1.0/dummy;
      }
    }
    stat->meansig=sum/((double)nval);
    stat->wmean=sumdatonsigsq/sum1onsigsq;
    stat->ewmean=1.0/sqrt(sum1onsigsq);
    if (efl!=NULL) stat->eflwmean=sqrt(sumeflonsigsq)/sum1onsigsq;
  } else {
    /** Use weights input by user **/
    for (i=0,sum=0.0; i<ndat; i++) {
      if (snull || sts[i]==1) {
	sum+=sig[i]; dummy=wgt[i]*wgt[i]; sumdatonsigsq+=dat[i]/dummy;
	sumsigonsigsq+=sig[i]*sig[i]/dummy/dummy;
	if (efl!=NULL) sumeflonsigsq+=efl[i]*efl[i]/dummy/dummy;
	sum1onsigsq+=1.0/dummy; 
      }
    }
    stat->meansig=sum/((double)nval);
    stat->wmean=sumdatonsigsq/sum1onsigsq;
    stat->ewmean=sqrt(sumsigonsigsq)/sum1onsigsq;
    if (efl!=NULL) stat->eflwmean=sqrt(sumeflonsigsq)/sum1onsigsq;
  }

  /* Calculate RMS of sigma-array, chisquared and reduced-chisquared
     around the weighted mean */
  for (i=0; i<ndat; i++) {
    if (snull || sts[i]==1) {
      dummy=sig[i]-stat->meansig; stat->rmssig+=dummy*dummy;
      /** Diverge here for use of expected fluctuations **/
      dummy=(efl==NULL) ? (dat[i]-stat->wmean)/sig[i] :
	(dat[i]-stat->wmean)/efl[i];
      stat->chisq+=dummy*dummy;
    }
  }
  stat->rmssig=(nval>1) ? sqrt(stat->rmssig/((double)(nval-1))) : 0.0;
  stat->rchisq=(nval>1) ? stat->chisq/((double)(nval-1)) : 0.0;
  
  /* Return here if medians are not required */
  if (opt==1) return 1;

  /* Find the median and semi-interquartile range of the data */
  if (opt==2) {
    if (!median(dat,sts,ndat,stat,0)) {
      nferrormsg("stats(): Error returned from median()\n\
\twhen taking median of dat array"); return 0;
    }
  } else if (opt==3) {
    if (!median(sig,sts,ndat,stat,0)) {
      nferrormsg("stats(): Error returned from median()\n\
\twhen taking median of sig array"); return 0;
    }
    stat->emed=stat->med; stat->esiqr=stat->siqr; stat->med=stat->siqr=0.0;
    if (!median(dat,sts,ndat,stat,0)) {
      nferrormsg("stats(): Error returned from median()\n\
\twhen taking median of dat array"); return 0;
    }
  }

  return 1;

}
