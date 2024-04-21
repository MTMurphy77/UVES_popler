/****************************************************************************
* Given a vaccum-heliocentric wavelength, calculate the value of the
* pixel index (double) which gives the corresponding
* air-wavelength. Function returns the double valued pixel index. This
* value is determined by knowledge of the vac-helio wavelengths
* corresponding to several (NINTP) positions along the pixel of
* interest (idx) and by using polynomial interpolation to derive more
* precise value.
*
* Opt=0 is for the wavelength scale of the object flux spectrum.
* Opt=1 is for the wavelength scale of the sky flux spectrum.
*
* Note: ranseed is the random number seed which determines the
* "offset" and "slope" parameters of the velocity distortions in
* UVES_wpol, when requested by the user. It should always come from
* the same spectrum (e.g. the first spectrum read in), but this can
* obviously be changed however the calling routine sees fit. Be sure
* that it's used consistently between calls of UVES_wpol and
* UVES_revwpol.
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "fit.h"
#include "memory.h"
#include "error.h"

double UVES_revwpol(spectrum *spec, int ord, int idx, double wl, long ranseed,
		    int opt, params *par) {

  double pix=0.0,dpix=0.0;
  double *x=NULL,*y=NULL;
  int    i=0;

  /* Allocate memory for interpolation data array */
  if ((x=darray(NINTP))==NULL)
    errormsg("UVES_revwpol(): Cannot allocate memory to x array of size %d",NINTP);
  if ((y=darray(NINTP))==NULL)
    errormsg("UVES_revwpol(): Cannot allocate memory to y array of size %d",NINTP);

  /* Fill x (i.e. wavelength) and y (i.e. pixel) interpolation arrays */
  if (opt) {
    x[0]=(idx) ? spec->ss.or[ord].vhrwl[idx-1] : spec->ss.or[ord].fvhlwl;
    x[NINTP-1]=spec->ss.or[ord].vhrwl[idx];
    if (wl<x[0] || wl>x[NINTP-1])
      warnmsg("UVES_revwpol(): Pixel %d of order %d in sky file\n\t%s\n\
\tdoes not include wavelength to interpolate, %lf",idx+1,ord+1,spec->ss.file,wl);
  } else {
    x[0]=(idx) ? spec->or[ord].vhrwl[idx-1] : spec->or[ord].fvhlwl;
    x[NINTP-1]=spec->or[ord].vhrwl[idx];
    if (wl<x[0] || wl>x[NINTP-1])
      warnmsg("UVES_revwpol(): Pixel %d of order %d in file\n\t%s\n\
\tdoes not include wavelength to interpolate, %lf",idx+1,ord+1,spec->file,wl);
  }
  y[0]=-0.5+(double)idx; y[NINTP-1]=y[0]+1.0;
  for (i=1; i<NINTP-1; i++) {
    y[i]=y[0]+(double)i/(double)(NINTP-1);
    x[i]=UVES_wpol(spec,ord,y[i],ranseed,opt,par);
  }

  /* Find pixel corresponding to input wavelength */
  if (!dpolint(x,y,NINTP,wl,&pix,&dpix))
    errormsg("UVES_revwpol(): Polynomial interpolation with %d points\n\
\tfailed when attempting to find pixel index corresponding to vac-helio\n\
\twavelength %lf",NINTP,wl);
  
  /* Clean up */
  free(x); free(y);

  return pix;

}
