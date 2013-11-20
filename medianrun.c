/****************************************************************************
MEDIANRUN: Function to calculate running median of a double array. The
median smoothing is calculated using a scale-number of pixels,
nfilt. The algorithm is taken from Paul Hewett's temcorr_4.f.

The program works by constructing an integer array ibuff which
contains the amplitude order of each pixel in dat rather than the
actual amplitude. An integer array is updated with the pixels within
the filter window and is then searched for the median value.  This
approach avoids the need to sort and is very fast. The sts input array
flags valid (sts=1) and invalid (sts!=1) pixels. If sts is passed as
NULL then all pixels in the dat array are treated as valid. The
resulting median array is passed back as med.

NOTE: nfilt must be an odd integer

****************************************************************************/

#include <stdlib.h>
#include "stats.h"
#include "sort.h"
#include "memory.h"
#include "error.h"

int medianrun(double *dat, double *med, int *sts, int ndat, int nfilt) {
 
  double  tmed=0.0;
  double  threemed[3];
  int     snull=0,nhalf=0,nval=0,oldval=0,newval=0,medlst1=0,medlst2=0;
  int     ne=0,ns=0,neo=0,nso=0,nmed=0,nmedhalf=0,nmede=0;
  int     i=0,j=0;
  int     *ibuff=NULL;
  dbleint *dibuff=NULL;

  /* Check nfilt is > 1 */
  if (nfilt<=1) {
    nferrormsg("medianrun(): nfilt (=%d) must be > 1",nfilt); return 0;
  }
  /* Check nfilt is odd */
  if (!(MOD(nfilt,2))) {
    nferrormsg("medianrun(): nfilt (=%d) must be odd",nfilt); return 0;
  }
  /* Check ndat>=nfilt */
  if (ndat<nfilt) {
    nferrormsg("medianrun(): Can not filter on a scale (nfilt=%d)\n\
\tlonger than the data array (ndat=%d)",nfilt,ndat); return 0;
  }

  /* Check if sts array was passed as null */
  if (sts==NULL) snull=1;

  /** Construct Running Median array **/
  if (nfilt==1) {
    /* The trivial case! */
    for (i=0; i<ndat; i++) med[i]=dat[i];
  } else if (nfilt==3) {
    /* Loop over pixels */
    for (i=0; i<ndat; i++) {
      /* Find start and end of median window */
      ns=MAX(0,(i-1)); ne=MIN((ndat-1),(i+1));
      /* Load valid data into threemed array while counting number of
	 valid pixels for us in median */
      for (j=ns,nmed=0; j<=ne; j++)
	if (snull || sts[j]==1) threemed[nmed++]=dat[j];
      /* Calculate median for those valid pixels */
      switch (nmed) {
      case 0: med[i]=0.0; break;
      case 1: med[i]=threemed[0]; break;
      case 2: med[i]=0.5*(threemed[0]+threemed[1]); break;
      case 3: 
	tmed=(MAX(threemed[0],threemed[1]));
	if (threemed[2]>=tmed) med[i]=tmed;
	else {
	  tmed=(MIN(threemed[0],threemed[1])); med[i]=(MAX(tmed,threemed[2]));
	}
	break;
      }
    }
  }
  else {
    /* For scales other than three pixels, define array ibuff that
       contains the order of each pixel value, rather than the actual
       value.  e.g. if the 10th pixel is the smallest value then
       ibuff[9]=1, if the 9th pixel is the next smallest then
       ibuff[8]=2, etc...*/
    /* Define half filter size */
    nhalf=nfilt/2;
    /* Allocate memory to arrays */
    if ((ibuff=iarray(ndat))==NULL)
      errormsg("medianrun(): Unable to allocate memory for ibuff array\n\
\tof size %d",ndat);
    if (!(dibuff=(dbleint *)malloc((size_t)((ndat)*sizeof(dbleint)))))
      errormsg("medianrun(): Unable to allocate memory for dibuff array\n\
\tof size %d",ndat);
    /* Fill temporary data and index arrays, leaving out invalid
       pixels and counting the number of valid ones */
    for (i=0,nval=0; i<ndat; i++) {
      if (snull || sts[i]==1) { dibuff[nval].a=dat[i]; dibuff[nval].i=i; nval++; }
    }
    /* Check to make sure some valid pixels were detected */
    if (!nval) { free(ibuff); free(dibuff); return 1; }
    /* Sort both double and integer parts of temporary array */
    qsort(dibuff,nval,sizeof(dbleint),qsort_dbleint);
    /* Fill in ibuff array and reset integer array part of dibuff */
    for (i=0; i<nval; i++) { ibuff[dibuff[i].i]=i+1; dibuff[i].i=0; }
    i=0; 
    /* Begin here when no prior valid median pixel has been defined
       (i.e. at the beginning of the data and after moving through any
       holes in the data equal to or larger than nfilt) */
    while (i<ndat) {
      /* Find first valid pixel and (re-)position starting and ending
	 pixels for filtering */
      while (i<ndat && !ibuff[i]) i++; if (i==ndat) break;
      i=(MAX(0,(i-nhalf))); ns=(MAX(0,(i-nhalf))); ne=(MIN((i+nhalf),(ndat-1)));
      /* Set integer array to 1 for all pixels in the filter window
	 centered on first pixel to be filtered */
      for (j=ns,nmed=0; j<=ne; j++)
	if (ibuff[j]) { dibuff[ibuff[j]-1].i=1; nmed++; }
      /* See if total number of available median points is odd or even */
      nmedhalf=nmed/2; nmede=(MOD((nmed+1),2));
      /* Find middle occupied value or two middle occupied values of
	 integer array to give filtered value for first pixel */
      if (!nmede) {
	medlst1=j=0; while (j<=nmedhalf) j+=dibuff[medlst1++].i;
	med[i]=dibuff[(--medlst1)].a;
      } else {
	medlst1=j=0; while (j<nmedhalf) j+=dibuff[medlst1++].i;
      	medlst2=medlst1--; while (j<=nmedhalf) j+=dibuff[medlst2++].i;
	med[i]=0.5*(dibuff[medlst1].a+dibuff[(--medlst2)].a);
      }
      /* Process remaining pixels. For each pixel delete the value
	 just moved out of the filter window from dibuff integer
	 array. Add the new value entering the filter to the integer
	 array. Then search integer array for the median value. The
	 direction of the search is specified by the size of the new
	 value relative to the last median value medlst and the value
	 just removed from the window */
      i++; while (i<ndat) {
	nso=ns; neo=ne; ns=(MAX(0,(i-nhalf))); ne=(MIN((i+nhalf),(ndat-1)));
	oldval=ibuff[nso]-1; newval=ibuff[ne]-1;
	if (!ibuff[nso] || ns==nso) {
	  if (ibuff[ne] && ne>neo) {
	    nmed++; nmede=(nmede) ? 0 : 1; dibuff[newval].i=1;
	    if (!nmede) {
	      if (newval>medlst2) medlst1=medlst2;
	      else if (newval>medlst1 && newval<medlst2) medlst1=newval;
	      med[i]=dibuff[medlst1].a;
	    } else {
	      if (newval>medlst1) {
		medlst2=medlst1+1; while (!dibuff[medlst2].i) medlst2++;
	      } else {
		medlst2=medlst1; medlst1--; while (!dibuff[medlst1].i) medlst1--;
	      }
	      med[i]=0.5*(dibuff[medlst1].a+dibuff[medlst2].a);
	    }
	  } else med[i]=med[i-1];
	} else {
	  if (ibuff[ne] && ne>neo) {
	    dibuff[oldval].i=0; dibuff[newval].i=1;
	    if (!nmede) {
	      if (oldval<=medlst1 && newval>medlst1) {
		medlst1++; while (!dibuff[medlst1].i) medlst1++;
		med[i]=dibuff[medlst1].a;
	      } else if (oldval>=medlst1 && newval<medlst1) {
		medlst1--; while (!dibuff[medlst1].i) medlst1--;
		med[i]=dibuff[medlst1].a;
	      } else med[i]=med[i-1];
	    } else {
	      if (oldval==medlst1 && newval<medlst2) {
		medlst1=medlst2-1; while (!dibuff[medlst1].i) medlst1--;
		med[i]=0.5*(dibuff[medlst1].a+dibuff[medlst2].a);
	      } else if (oldval==medlst2 && newval>medlst1) {
		medlst2=medlst1+1; while (!dibuff[medlst2].i) medlst2++;
		med[i]=0.5*(dibuff[medlst1].a+dibuff[medlst2].a);
	      } else if (oldval<medlst1 && newval>medlst1 && newval<medlst2) {
		medlst1=newval; med[i]=0.5*(dibuff[medlst1].a+dibuff[medlst2].a);
	      } else if (oldval>medlst2 && newval>medlst1 && newval<medlst2) {
		medlst2=newval; med[i]=0.5*(dibuff[medlst1].a+dibuff[medlst2].a);
	      } else if (oldval<=medlst1 && newval>medlst2) {
		medlst1=medlst2++; while (!dibuff[medlst2].i) medlst2++;
		med[i]=0.5*(dibuff[medlst1].a+dibuff[medlst2].a);
	      } else if (oldval>=medlst2 && newval<medlst1) {
		medlst2=medlst1--; while (!dibuff[medlst1].i) medlst1--;
		med[i]=0.5*(dibuff[medlst1].a+dibuff[medlst2].a);
	      } else med[i]=med[i-1];
	    }
	  } else {
	    if (!(--nmed)) break; nmede=(nmede) ? 0 : 1; dibuff[oldval].i=0;
	    if (!nmede) {
	      if (oldval<=medlst1) medlst1=medlst2;
	      med[i]=dibuff[medlst1].a;
	    } else {
	      if (oldval>medlst1) {
		medlst2=medlst1--; while (!dibuff[medlst1].i) medlst1--;
	      } else if (oldval<medlst1) {
		medlst2=medlst1+1; while (!dibuff[medlst2].i) medlst2++;
	      } else {
		medlst2=medlst1--; while (!dibuff[medlst1].i) medlst1--;
		medlst2++; while (!dibuff[medlst2].i) medlst2++;
	      }
	      med[i]=0.5*(dibuff[medlst1].a+dibuff[medlst2].a);
	    }
	  }
	}
	/* Advance to next pixel */
	i++;
      }
    }
    /* Clean up */
    free(ibuff); free(dibuff);
  }

  return 1;

}
