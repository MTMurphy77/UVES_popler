/****************************************************************************
* Set the wavelength scale of the combined spectrum
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int UVES_set_wavelen_scale(spectrum *spec, int nspec, cspectrum *cspec,
			   params *par) {

  double cwls=INFIN,cwle=0.0,cdisp=INFIN; /* Combined spectrum wavelength params */
  double v_saw=0.0,dv=0.0;
  int    i=0,j=0,k=0,l=0;
  int    v_distort=0;
  int    *ibuf=NULL;
  long   ranseed=0;

  /* Introduce a +/- v_saw m/s saw-tooth wavelength distortion
     to each order like that found by Griest et al. 2010 */
  /* Set v_distort to 1 to apply distortions, 0 to NOT apply distortions. */
  v_distort=0;
  // For Keck
  v_saw=200.0;
  // For VLT
  // v_saw=100.0;

  /* Just for convenience and formatting below, copy random number
     seed from first spectrum read in for applying velocity
     distortions in UVES_wpol calls (nothing to do with velocity
     distortions discussed above) */
  ranseed=spec[0].distort_seed;
  
  /* Use wavelength calibration polynomial to set the vaccuum-heliocentric
     scale for each order, determine the starting and ending wavelengths of
     combined spectrum and determine the dispersion of the combined
     spectrum */
  for (i=0; i<nspec; i++) {
    /* Loop over orders in each spectrum */
    for (j=0; j<spec[i].nor; j++) {
      /* Only operate on useful orders */
      if (spec[i].or[j].nuse>=MINUSE) {
	/* Test dispersion to see if it is reversed or not */
	spec[i].or[j].revdisp=0;
	spec[i].or[j].vhwl[0]=UVES_wpol(&(spec[i]),j,0.0,ranseed,par);
	spec[i].or[j].vhrwl[0]=UVES_wpol(&(spec[i]),j,0.5,ranseed,par);
	if (spec[i].or[j].vhrwl[0]<spec[i].or[j].vhwl[0]) {
	  /* Dispersion must be reversed */
	  spec[i].or[j].revdisp=1;
	  /* When dispersion is reversed input file, first rearrange
	     the flux, error and status arrays */
	  /* Allocate memory for temporary status array */
	  if ((ibuf=iarray(spec[i].or[j].np))==NULL)
	    errormsg("UVES_set_wavelen_scale(): Cannot allocate memory for\n\
\tibuf array of size %d",spec[i].or[j].np);
	  for (k=0,l=spec[i].or[j].np-1; k<spec[i].or[j].np; k++,l--) {
	    spec[i].or[j].vhwl[k]=spec[i].or[j].fl[l];
	    spec[i].or[j].vhrwl[k]=spec[i].or[j].er[l];
	    ibuf[k]=spec[i].or[j].st[l];
	  }
	  for (k=0; k<spec[i].or[j].np; k++) {
	    spec[i].or[j].fl[k]=spec[i].or[j].vhwl[k];
	    spec[i].or[j].er[k]=spec[i].or[j].vhrwl[k];
	    spec[i].or[j].st[k]=ibuf[k];
	  }
	  /* Clean up */
	  free(ibuf);
	  /* Reinstate wavelength information for first pixel */
	  spec[i].or[j].vhwl[0]=UVES_wpol(&(spec[i]),j,0.0,ranseed,par);
	  spec[i].or[j].vhrwl[0]=UVES_wpol(&(spec[i]),j,0.5,ranseed,par);
	}
	/* Set left-hand vac. helio. wavelength of first pixel */
	spec[i].or[j].fvhlwl=UVES_wpol(&(spec[i]),j,-0.5,ranseed,par);
	/* For CPL UVES spectra, convert the resolution array from
	   pixel space into velocity space */
	if (spec->ftype==FTUVES && spec->fvers==1)
	  spec[i].or[j].res[0]*=C_C_K*(spec[i].or[j].vhrwl[0]-spec[i].or[j].fvhlwl)/
	    spec[i].or[j].vhwl[0];
	/* Set the wavelengths of the middles and right-hand edges of the pixels */
	for (k=1; k<spec[i].or[j].np; k++) {
	  /* Following line was commented out for version 0.64 (20th
	     June 2013). It's possible that "FTIESI" was supposed to
	     be "FTHIRX". In any case, we now read in the raw
	     wavelength array for MAGE and HIRES REDUX files and use
	     them in UVES_wpol() (instead of the mistaken use of the
	     vac. helio. wav. arrays) to define the final wavelength
	     scale. So the following line should not be necessary. */
	  // if (spec->ftype!=FTMAGE && spec->ftype!=FTIESI)
	  spec[i].or[j].vhwl[k]=UVES_wpol(&(spec[i]),j,(double)k,ranseed,par);
	  spec[i].or[j].vhrwl[k]=UVES_wpol(&(spec[i]),j,0.5+(double)k,ranseed,par);
	  if (par->linear)
	    cdisp=(MIN(cdisp,(spec[i].or[j].vhwl[k]-spec[i].or[j].vhwl[k-1])));
	  else
	    cdisp=(MIN(cdisp,(log10(spec[i].or[j].vhwl[k]/spec[i].or[j].vhwl[k-1]))));
	  /* For CPL UVES spectra, convert the resolution array from
	     pixel space into velocity space */
	  if (spec->ftype==FTUVES && spec->fvers==1)
	    spec[i].or[j].res[k]*=C_C_K*
	      (spec[i].or[j].vhrwl[k]-spec[i].or[j].vhrwl[k-1])/spec[i].or[j].vhwl[k];
	}
	/* Introduce a +/- v_saw m/s saw-tooth wavelength distortion
	   to each order like that found by Griest et al. 2010 */
	if (v_distort) {
	  for (k=0; k<spec[i].or[j].np; k++) {
	    dv=1.0+(v_saw-4.0*v_saw*fabs(0.5-(double)k/(double)spec[i].or[j].np))/C_C;
	    spec[i].or[j].vhwl[k]*=dv; spec[i].or[j].vhrwl[k]*=dv; 
	  }
	}
	/* Calculate number of valid pixels */
	k=0; while (spec[i].or[j].st[k]<1) k++;
	cwls=(k) ? (MIN(cwls,spec[i].or[j].vhrwl[k-1])) :
	  (MIN(cwls,spec[i].or[j].fvhlwl));
	k=spec[i].or[j].np-1; while (spec[i].or[j].st[k]<1) k--;
	cwle=MAX(cwle,spec[i].or[j].vhrwl[k]);
      }
    }
  }
  if (cwls==INFIN || cwle==0.0)
    errormsg("UVES_set_wavelen_scale(): No valid pixels found");

  /* Set dispersion */
  if (par->linear) {
    cspec->dv=0.0;
    cspec->dwl=(par->disp<DRNDTOL) ? (par->disp=cdisp) : (cdisp=par->disp);
  } else {
    cspec->dwl=0.0;
    cspec->dv=(par->disp<DRNDTOL) ? cdisp : (cdisp=log10(par->disp/C_C_K+1.0));
    if (par->disp<DRNDTOL) par->disp=C_C_K*(pow(10.0,cspec->dv)-1.0);
  }

  /* Determine number of pixels in final combined spectrum */
  if (par->linear) cspec->np=(int)((cwle-cwls)/cdisp)+1;
  else cspec->np=(int)(log10(cwle/cwls)/cdisp)+1;

  /* Fill first pixel's left edge wavelength */
  cspec->flwl=cwls;
  /* TEMPORARY - For introducing an offset in the binning used */
  // cspec->flwl*=pow(10.0,-0.50*cspec->dv);

  /* Initialize the combined spectrum and set its wavelength scale */
  if (!UVES_init_cspec(cspec,par,0)) {
    nferrormsg("UVES_set_wavelen_scale(): Error returned from\n\
\tUVES_init_cspec()"); return 0;
  }

  return 1;

}
