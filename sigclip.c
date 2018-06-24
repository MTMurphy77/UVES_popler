/****************************************************************************
SIGCLIP: Function to sigma-clip a given array of data

The clip array is an array of 0's or 1's showing whether each element
of dat array has been clipped or not respectively. The statistics for
the remaining array of points are returned in stat. Input elements in
the clip array which are not equal to 1 will cause the corresponding
data elements to be ignored. Other elements in the clip array are
intialized to 1.

Opt=0 will use the statistical uncertainty (err) array as the "sigma"
in the clipping algorithm (but not the weighted mean calculation),
while any other value for opt will cause the expected fluctuation
(efl) array to be used.

****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "stats.h"
#include "memory.h"
#include "error.h"

int sigclip(double *dat, double *err, double *efl, double *wgt, int ndat,
	    double clipsig, statset *stat, int *clip, int opt) {

  double   sig=0.0,sign=0.0;
  int      sigidx=0,nval=0;
  int      i=0,j=0;

  /* Initialize clip array */
  for (i=0,j=0,nval=0; i<ndat; i++) if (clip[i]==1) { nval++; j=i; }

  /* Return if bad array length passed */
  if (nval<1) {
    nferrormsg("sigclip(): Only %d data points passed have clip=1",nval);
    return 0;
  }
 
  /* Deal with simple case */
  if (nval==1) {
    if (!stats(&(dat[j]),&(err[j]),&(efl[j]),&(wgt[j]),&(clip[j]),nval,1,stat)) {
      nferrormsg("sigclip(): Error returned from stats()\n\
\twhen calculating stats for %d valid data value"); return 0;
    }
    return 1;
  }

  /* Iterate, removing most deviant points 'til none are clipped from
     weighted mean */
  while (nval>1) {
    /* Calculate weighted mean */
    if (!stats(dat,err,efl,wgt,clip,ndat,1,stat)) {
      nferrormsg("sigclip(): Error returned from stats()\n\
\twhen calculating stats for %d valid data values",nval); return 0;
    }
    /* Find most deviant point */
    i=0; while (clip[i]!=1) i++;
    if (!opt) {
      sig=fabs(dat[i]-stat->wmean)/err[i];
      sigidx=i; for (j=i+1; j<ndat; j++) {
	if (clip[j]==1 && (sign=fabs(dat[j]-stat->wmean)/err[j])>sig) {
	  sig=sign; sigidx=j;
	}
      }
    } else {
      sig=fabs(dat[i]-stat->wmean)/efl[i];
      sigidx=i; for (j=i+1; j<ndat; j++) {
	if (clip[j]==1 && (sign=fabs(dat[j]-stat->wmean)/efl[j])>sig) {
	  sig=sign; sigidx=j;
	}
      }
    }
    /* Clip out most deviant point */
    if (sig>clipsig) { clip[sigidx]=0; nval--; }
    else break;
  }

  /* Do final statistics */
  if (!stats(dat,err,efl,wgt,clip,ndat,1,stat)) {
    nferrormsg("sigclip(): Error returned from stats()\n\
\twhen calculating stats for %d valid dafta values",nval); return 0;
  }

  return 1;

}
