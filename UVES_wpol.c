/****************************************************************************
* Calculate the vaccuum-heliocentric wavelength of a given pixel from the
* wavelength calibration polynomial.
*
* Note: ranseed is the random number seed which determines the
* "offset" and "slope" parameters of the velocity distortions, when
* requested by the user. It should always come from the same spectrum
* (e.g. the first spectrum read in), but this can obviously be changed
* however the calling routine sees fit. The random number seed
* controlling the "intra-order distortion" is the one attached to the
* structure of the spectrum whose wavelength scale is being determined
* in this routine.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "fit.h"
#include "stats.h"
#include "memory.h"
#include "const.h"
#include "error.h"

double UVES_wpol(spectrum *spec, int ord, double idx, long ranseed, params *par) {

  double airwl=0.0,vacwl=0.0,didx=0.0,ndidx=0.0,dndidx=0.0;
  double dwl=0.0;
  /* Velocity distortion parameter limits in m/s. The range these
     parameters take after (uniform) randomization is +/- the values
     set here. The velocity offset is applied at 5000 Angstrom and the
     slope here is defined in units of (m/s)/Ang. over the range
     3000-10000Ang. */
  double vd_offset=50.0,vd_slope=0.02,vd_echamp=50.0;
  double vd_x=0.0,vd=0.0;
  double *x=NULL,*y=NULL;
  int    sidx=0;
  int    i=0,j=0;
  long   vd_idum=0;

  /* Wavelength calibration polynomials are usually 1-referenced,
     rather than 0-referenced (as in these C-codes). We therefore need
     to add 1 to the input index here, where appropriate */
  if (spec->ftype==FTUVES && spec->fvers==1) didx=idx+1.0-spec->or[ord].wpol[NWPOL-1];
  else if (spec->ftype==FTIRLS || spec->ftype==FTHIRX || spec->ftype==FTIESI) didx=idx;
  else if (spec->ftype==FTMAKE) didx=idx+spec->or[ord].idxoff;
  else didx=idx+1.0;

  /* If the dispersion (as derived from the wavelength polynomial) of
     the echelle order of interest is reversed then we need to
     transform the input index appropriately */
  if (spec->or[ord].revdisp) didx=(double)(spec->or[ord].np)-didx+1.0;

  /* Use the wavelength polynomial coefficients to calculate air-wavelength */
  switch (spec->ftype) {
  case FTUVES:
    for (i=spec->or[ord].nwpol-1; i>=0; i--) airwl=spec->or[ord].wpol[i]+didx*airwl;
    if (!spec->fvers) airwl*=1.e4;
    else if (spec->fvers==1) airwl+=spec->or[ord].wpol[NWPOL-2];
    break;
  case FTIRAF:
    ndidx=(2.0*didx-(double)(spec->or[ord].np+1))/(double)(spec->or[ord].np-1);
    dndidx=2.0/(double)(spec->or[ord].np-1);
    if (ndidx<-1.0 || ndidx>1.0) {
      /* Allocate memory for interpolation data array */
      if ((x=darray(NINTP))==NULL)
	errormsg("UVES_wpol(): Cannot allocate memory to x array of size %d",NINTP);
      if ((y=darray(NINTP))==NULL)
	errormsg("UVES_wpol(): Cannot allocate memory to y array of size %d",NINTP);
      if (ndidx<-1.0) {
	for (i=0; i<NINTP; i++) {
	  x[i]=-1.0+dndidx*(double)i/(double)(NINTP-1);
	  y[i]=chebyshev_eval(x[i],spec->or[ord].wpol,spec->or[ord].nwpol);
	}
      } else {
	for (i=0; i<NINTP; i++) {
	  x[i]=1.0-dndidx*(1.0-(double)i/(double)(NINTP-1));
	  y[i]=chebyshev_eval(x[i],spec->or[ord].wpol,spec->or[ord].nwpol);
	}
      }
      /* Find pixel corresponding to input wavelength */
      if (!dpolint(x,y,NINTP,ndidx,&airwl,&dwl))
	errormsg("UVES_wpol(): Polynomial interpolation with %d points\n\
\tfailed when attempting to find wavelength corresponding to pixel %lf\n\
\tin order %d of file %d,\n\t%s",NINTP,idx,ord+1,spec->id,spec->file);
      /* Clean up */
      free(x); free(y);
    } else airwl=chebyshev_eval(ndidx,spec->or[ord].wpol,spec->or[ord].nwpol);
    break;
  case FTMAKE:
    for (i=spec->or[ord].nwpol-1; i>=0; i--) airwl=spec->or[ord].wpol[i]+didx*airwl;
    break;
  case FTIRLS:
    /* Linear wavelength solution stored in wpol parameters */
    /* Dispersion: w[i]=CRVAL1+([i+1]-CRPIX1)*CD1_1 */
    for (i=spec->or[ord].nwpol-1; i>=0; i--) airwl=spec->or[ord].wpol[i]+didx*airwl;
    break;
  case FTHIRX:
    if (par->thar<=1) {
      /* Allocate memory for interpolation data array */
      if ((x=darray(NINTP))==NULL)
	errormsg("UVES_wpol(): Cannot allocate memory to x array of size %d",NINTP);
      if ((y=darray(NINTP))==NULL)
	errormsg("UVES_wpol(): Cannot allocate memory to y array of size %d",NINTP);
      sidx=(MAX(0,(NINT((MIN((didx-(double)NINTP/2.0),
			     ((double)(spec->or[ord].np-NINTP-1))))))));
      for (i=0,j=sidx; i<NINTP; i++,j++) { x[i]=(double)j; y[i]=spec->or[ord].wl[j]; }
      /* Find wavelength corresponding to input pixel */
      if (!dpolint(x,y,NINTP,didx,&airwl,&dwl))
	errormsg("UVES_wpol(): Polynomial interpolation with %d points\n\
\tfailed when attempting to find wavelength corresponding to pixel %lf\n \
\tin order %d of file %d,\n\t%s",NINTP,idx,ord+1,spec->id,spec->file);
      /* Clean up */
      free(x); free(y);
    } else if (par->thar==2) {
      /* Stuff below is for using the wavelength polynomial read in from
	 the ThAr files themselves. See comment in UVES_r2Dspec_hirx() */
      ndidx=2.0*(didx-spec->or[ord].wpol[NWPOL-2])/spec->or[ord].wpol[NWPOL-1];
      airwl=legendre_eval(ndidx,spec->or[ord].wpol,spec->or[ord].nwpol);
    }
    break;
  case FTESOM:
    /* Linear wavelength solution stored in wpol parameters */
    /* Dispersion: w[i]=CRVAL1+([i+1]-CRPIX1)*CDELT1 */
    for (i=spec->or[ord].nwpol-1; i>=0; i--) airwl=spec->or[ord].wpol[i]+didx*airwl;
    break;
  case FTMAGE:
    /* Wavelength scale is provided as a look-up table. The scale is
       very close to being a log-linear scale, but is not quite
       log-linear. To determine wavelength values between pixels, use
       local log-linear dispersion, log10(w[i+1]/w[i]) */
    i=(MIN((MAX((int)idx,0)),(spec->or[ord].np-2)));
    dwl=log10(spec->or[ord].wl[i+1]/spec->or[ord].wl[i]);
    airwl=spec->or[ord].wl[i]*pow(10.0,(idx-(double)i)*dwl);
    break;
  case FTIESI:
    /* Linear wavelength solution stored in wpol parameters */
    /* Dispersion: w[i]=wstart_ord+i*dispersion_ord */
    for (i=spec->or[ord].nwpol-1; i>=0; i--) airwl=spec->or[ord].wpol[i]+didx*airwl;
    break;
  }

  /* If user wants to apply a distortion, with randomized parameters,
     to the wavelength scale, do it here. The velocity shift applied
     to a given wavelength is (i) linearly dependent on wavelength,
     with random offset at 500nm and random slope, with the random
     seed coming from the first spectrum read in, (ii) quadratically
     dependent on pixel position within the echelle order, with a
     random (uniform about zero; from that spectrum's random seed)
     velocity amplitude of the distortion from the middle of the order
     to its edge.
  */
  if (par->distort) {
    vd_idum=ranseed;
    vd_offset*=2.0*(ran(&vd_idum)-0.5); vd_slope*=2.0*(ran(&vd_idum)-0.5);
    vd_idum=spec->distort_seed;
    vd_echamp*=2.0*(ran(&vd_idum)-0.5);
    vd_x=didx/(double)spec->or[ord].np-0.5;
    vd=vd_offset+vd_slope*(airwl-5000.0)+vd_echamp*(1.0-8.0*vd_x*vd_x);
    airwl+=airwl*vd/C_C;
  }

  /* Don't do anything else if only using ThAr frames. Also,
     HIRES_REDUX files are already in vacuum heliocentric frame, as
     are MAGE files; IRAF-ESI files are in vacuum, but not
     heliocentric frame. Just ensure that any user-supplied velocity
     shifts are applied. */
  if (par->thar==2 || spec->ftype==FTHIRX || spec->ftype==FTMAGE) {
    vacwl=airwl;
    // vacwl+=vacwl*spec->vshift/C_C_K;
    vacwl+=vacwl*(spec->vshift+(vacwl-spec->refwav)*spec->vslope/1000.0)/C_C_K;
    return vacwl;
  }

  /* If an air-vac wavelength conversion is required, do it now */
  if (!par->vacwl && spec->ftype!=FTIESI) {
    /* Solve Edlen formula iteratively to find the vacuum wavelength */
    if ((vacwl=edlen_a2v(airwl))==0.0)
      errormsg("UVES_wpol(): Problem converting air-wavlength %lf\n\
\tto vacuum wavelength scale in file\n\t%s",airwl,spec->file);
  } else vacwl=airwl;

  /* If conversion to the heliocentric frame is required, do it now */
  if (!par->helio) {
    /* Then convert to heliocentric frame */
    vacwl+=vacwl*spec->vhel/C_C_K;
  }

  /* Apply user-supplied velocity shift for this spectrum */
  // vacwl+=vacwl*spec->vshift/C_C_K;
  vacwl+=vacwl*(spec->vshift+(vacwl-spec->refwav)*spec->vslope/1000.0)/C_C_K;

  return vacwl;

}
