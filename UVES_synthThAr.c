/****************************************************************************
* Generate a synthetic emission line (or ThAr) calibration spectrum as
* a series of equal height, equal width (in velocity space)
* Gaussians. The existing ThAr spectra are overwritten in this
* process. No noise is put on the new spectra but the error array is
* reset to a small value.
****************************************************************************/

#include <math.h>
#include "UVES_popler.h"
#include "gamm.h"
#include "const.h"
#include "error.h"

int UVES_synthThAr(spectrum *spec, params *par) {

  double   wlg=0.0,pref=0.0,sigl=0.0,sigr=0.0;
  /* double   wlg2=0.0,pref2=0.0,sigl2=0.0,sigr2=0.0; */
  int      fval=0,lval=0; /* First and last valid pixels in raw order */
  int      i=0,l=0;

  /* Loop over all orders */
  for (l=0; l<spec->nor; l++) {
    /* Check first to see if order is useful */
    if (spec->or[l].nuse>=MINUSE) {
      /* Find first valid pixel in raw array */
      i=0; while (spec->or[l].st[i]<1) i++; fval=i;
      /* Find last valid pixel in raw array */
      i=spec->or[l].np-1; while (spec->or[l].st[i]<1) i--; lval=i;
      /* Loop over valid pixels */
      for (i=fval; i<=lval; i++) {
	/* Is this pixel near enough to a Gaussian for us to calculate the flux? */
	wlg=2.0*(double)((int)(spec->or[l].vhwl[i]/2.0));
	wlg=(spec->or[l].vhwl[i]-wlg>1.0) ? wlg+2.0 : wlg; 
	/* Initialise flux and error */
	spec->or[l].fl[i]=100.0; spec->or[l].er[i]=0.1;
	if (C_C_K*fabs(spec->or[l].vhwl[i]-wlg)/wlg<20.0) {
	  pref=70000.0*C_FWHMSIG/spec->or[l].vhwl[i]/C_SQRT2;
	  sigl=(i) ? pref*(spec->or[l].vhrwl[i-1]-wlg) :
	    pref*(spec->or[l].fvhlwl-wlg);
	  sigr=pref*(spec->or[l].vhrwl[i]-wlg);
	  spec->or[l].fl[i]+=10000.0*(erffn(sigr)-erffn(sigl));
	  spec->or[l].er[i]+=1.e-2*sqrt(spec->or[l].fl[i]-100.0);
	  /*
	  wlg2=wlg*(1.0-0.4/C_C_K);
	  pref2=110000.0*C_FWHMSIG/spec->or[l].vhwl[i]/C_SQRT2;
	  sigl2=(i) ? pref2*(spec->or[l].vhrwl[i-1]-wlg2) :
	    pref2*(spec->or[l].fvhlwl-wlg2);
	  sigr2=pref2*(spec->or[l].vhrwl[i]-wlg2);
	  spec->or[l].fl[i]+=10000.0*(erffn(sigr2)-erffn(sigl2));
	  */
	}
      }
    }
  }
      
  return 1;

}
