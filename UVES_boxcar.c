/****************************************************************************
* Boxcar-smooth a spectrum within a sliding window of nwid pixels or a
* velocity vwidth of vwid km/s. Option 1 determines whether the
* weighted mean flux is taken within the sliding window or just the
* mean (1 & 0 respectively). Option 2 determines whether nwid or vwid
* is used (0 & 1 respectively).
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "stats.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int UVES_boxcar(double *wl, double *fl, double *er, int ndat, int nwid,
		double vwid, double *mean, double *rms, int opt1, int opt2) {

  double   disp=0.0;
  double   *dat=NULL,*err=NULL;
  int      nsm=1,sidx=0,nval=0;
  int      i=0,j=0,k=0;
  statset  stat;

  /* Make sure it's reasonable to smooth the given array */
  if (ndat<2) {
    nferrormsg("UVES_boxcar(): Cannot smooth data array only\n\
\t%d pixels long",ndat);
    return 0;
  }

  if (!opt2) nsm=nwid;

  /* Allocate initial memory for statistics arrays */
  if ((dat=darray(nsm))==NULL)
    errormsg("UVES_boxcar(): Cannot allocate memory for\n\
\tdat array of length %d",nsm);
  if ((err=darray(nsm))==NULL)
    errormsg("UVES_boxcar(): Cannot allocate memory for\n\
\terr array of length %d",nsm);

  /* Loop over pixels */
  for (i=0; i<ndat; i++) {
    if (opt2) {
      disp=(i==ndat-1) ? C_C_K*(wl[ndat-2]-wl[ndat-1])/wl[ndat-1] :
	C_C_K*(wl[i+1]-wl[i])/wl[i];
      nsm=1+2*(int)(0.5*vwid/disp);
      /* Reallocate memory for statistics arrays */
      if (!(dat=(double *)realloc(dat,(size_t)(nsm*sizeof(double)))))
	errormsg("UVES_boxcar(): Cannot reallocate memory for\n\
\tdat array of length %d",nsm);
      if (!(err=(double *)realloc(err,(size_t)(nsm*sizeof(double)))))
	errormsg("UVES_boxcar(): Cannot reallocate memory for\n\
\terr array of length %d",nsm);
    }
    sidx=i-nsm/2; sidx=(MIN((MAX(0,sidx)),(ndat-nsm)));
    for (j=sidx,k=0,nval=0; k<nsm; j++,k++)
      if (er[j]>0.0) { dat[nval]=fl[j]; err[nval++]=er[j]; }
    if (!nval) { mean[i]=0.0; rms[i]=0.0; }
    else {
      if (!stats(dat,err,NULL,NULL,NULL,nval,opt1,&stat))
	errormsg("UVES_boxcar(): Error returned from stats()");
      mean[i]=(opt1) ? stat.wmean : stat.mean; rms[i]=stat.rms;
    }
  }

  /* Clean up */
  free(dat); free(err);

  return 1;

}
