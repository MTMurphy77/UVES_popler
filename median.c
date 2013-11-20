/****************************************************************************
MEDIAN: Function to calculate median of a double array

There are usually extensive considerations as to how one might
implement a median-finding routine. The following routine bypasses all
that and just copies the array to a temporary array, sorts it and then
finds the middle value. Simple. Other methods will save time and will
be more convenient. This routine is just the simplest modular
solution. I also calculate the 68% semi-interquartile range (i.e. half
the range of data around the median which contains 68% of the
values). The status array, sts, should be 1 for valid pixels and any
other value for invalid pixels. If the status array is passed as NULL
then all pixels are assumed to be valid.

opt = 0 : If ndat is even then the two middle values are averaged to
          find the median.
opt = 1 : If ndat is even then the lower of the two middle values is
          returned as the median.
opt = 2 : If ndat is even then the higher of the two middle values is
          returned as the median.

****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "sort.h"
#include "stats.h"
#include "memory.h"
#include "error.h"

int median(double *dat, int *sts, int ndat, statset *stat, int opt) {

  double       *dbuf=NULL;
  int          snull=0,nval=0,sidx=0;
  register int i=0;

  /* Determine iff sts was passed as null or not */
  if (sts==NULL) snull=1;

  /* Return with non-fatal error if number of elements zero */
  if (ndat==0) {
    nferrormsg("median(): Number of data points passed is zero"); return 0;
  } else if (ndat==1) {
    stat->med=(snull || sts[0]==1) ? dat[0] : 0.0; stat->siqr=0.0; return 1;
  }
  if (ndat==2) {
    if (snull || (sts[0]==1 && sts[1]==1)) {
      if (opt==0) stat->med=0.5*(dat[0]+dat[1]);
      else if (opt==1) stat->med=(dat[0]<dat[1]) ? dat[0] : dat[1];
      else if (opt==2) stat->med=(dat[0]>dat[1]) ? dat[0] : dat[1];
      stat->siqr=0.5*fabs(dat[0]-dat[1]); return 1;
    } else if (sts[0]==1) {
      stat->med=dat[0]; stat->siqr=0.0; return 1;
    } else if (sts[1]==1) {
      stat->med=dat[1]; stat->siqr=0.0; return 1;
    } else { stat->med=stat->siqr=0.0; return 1; }
  }

  /* Allocate memory for temporary array */
  if ((dbuf=darray(ndat))==NULL) {
    nferrormsg("median(): Cannot allocate memory for temporary\n\
\tarray of size %d",ndat);
    return 0;
  }

  /* Copy and sort the array */
  for (i=0,nval; i<ndat; i++) if (snull || sts[i]==1) dbuf[nval++]=dat[i];
  if (nval) {
    /* Sort the array */
    qsort(dbuf,nval,sizeof(double),qsort_darray);
    /* Find median value and semi-interquartile range */
    sidx=nval/2; i=(int)(MED_QR*(double)nval);
  }

  /* Find the median */
  if (!nval) stat->med=stat->siqr=0.0;
  else if (nval==1) { stat->med=dbuf[0]; stat->siqr=0.0; }
  else if (nval==2) {
    if (opt==0) stat->med=0.5*(dbuf[0]+dbuf[1]);
    else if (opt==1) stat->med=(dbuf[0]<dbuf[1]) ? dbuf[0] : dbuf[1];
    else if (opt==2) stat->med=(dbuf[0]>dbuf[1]) ? dbuf[0] : dbuf[1];
    stat->siqr=0.5*fabs(dbuf[0]-dbuf[1]);
  } else if (isodd(nval)) {
    stat->med=dbuf[sidx]; stat->siqr=0.5*(dbuf[sidx+i]-dbuf[sidx-i]);
  } else {
    if (opt==0) stat->med=0.5*(dbuf[sidx-1]+dbuf[sidx]);
    else if (opt==1) stat->med=dbuf[sidx-1];
    else if (opt==2) stat->med=dbuf[sidx];
    stat->siqr=0.25*(-dbuf[sidx-i-1]-dbuf[sidx-i]+dbuf[sidx+i-1]+dbuf[sidx+i]);
  }

  /* Clean up */
  free(dbuf);

  return 1;

}
