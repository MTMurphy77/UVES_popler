/****************************************************************************
* For each order in a ThAr spectrum, create an error array by sliding
* a window across the order and determining the minimum RMS over the
* order and the mean flux at the position where the RMS is
* minimum. The error for each pixel is then just the minimum RMS times
* the square-root of the ratio of the flux in the pixel and the mean
* flux where the RMS is minimum.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

int UVES_thar_sigarray(spectrum *spec, params *par) {

  double   minrms=0.0,minmean=0.0;
  double   *mean=NULL;
  int      minidx=0;
  int      i=0,j=0;

  /* Loop over all orders */
  for (i=0; i<spec->nor; i++) {
    /* Check first to see if order is useful */
    if (spec->or[i].nuse>=MINUSE) {
      /* Allocate memory for temporary mean array */
      if ((mean=darray(spec->or[i].np))==NULL)
	errormsg("UVES_thar_sigarray(): Cannot allocate memory for\n\
\tmean array of length %d",spec->or[i].np);
      /* Pass sliding window along order to determine the RMS for
	 individual pixels */
      if (par->thar==1) {
	/* For the case when the ThAr spectra are held in their own
	   dedicated arrays */
	if (!UVES_boxcar(NULL,spec->or[i].th,spec->or[i].th,spec->or[i].np,
			 par->nordsig,0.0,mean,spec->or[i].ter,0,0)) {
	  nferrormsg("UVES_thar_sigarray(): Error returned from UVES_boxcar()\n\
\twhen operating on order %d",i+1); return 0;
	}
	/* Find minimum RMS across order and mean at that position */
	j=0;
	while (mean[j]<=DRNDTOL || spec->or[i].ter[j]<=DRNDTOL) j++;
	minidx=j; minrms=spec->or[i].ter[j];
	for (j=minidx+1; j<spec->or[i].np; j++)
	  if (mean[j]>DRNDTOL && spec->or[i].ter[j]>DRNDTOL &&
	      spec->or[i].ter[j]<minrms) { minrms=spec->or[i].ter[j]; minidx=j; }
	minmean=(MAX(mean[minidx],minrms));
	/* Model for the noise in every pixel */
	for (j=0; j<spec->or[i].np; j++)
	  spec->or[i].ter[j]=minrms*sqrt((MAX(spec->or[i].th[j],minmean))/
					 minmean);
      } else {
	/* For the case when the ThAr spectra are being combined and are
	   held in the object arrays */
	if (!UVES_boxcar(NULL,spec->or[i].fl,spec->or[i].fl,spec->or[i].np,
			 par->nordsig,0.0,mean,spec->or[i].er,0,0)) {
	  nferrormsg("UVES_thar_sigarray(): Error returned from UVES_boxcar()\n\
\twhen operating on order %d",i+1); return 0;
	}
	/* Find minimum RMS across order and mean at that position */
	j=0;
	while (mean[j]<=DRNDTOL || spec->or[i].er[j]<=DRNDTOL) j++;
	minidx=j; minrms=spec->or[i].er[j];
	for (j=minidx+1; j<spec->or[i].np; j++)
	  if (mean[j]>DRNDTOL && spec->or[i].er[j]>DRNDTOL &&
	      spec->or[i].er[j]<minrms) { minrms=spec->or[i].er[j]; minidx=j; }
	minmean=(MAX(mean[minidx],minrms));
	/* Model for the noise in every pixel */
	for (j=0; j<spec->or[i].np; j++)
	  spec->or[i].er[j]=minrms*sqrt((MAX(spec->or[i].fl[j],minmean))/minmean);
      }
      /* Clea up */
      free(mean);
    }
  }
  
  return 1;

}
