/****************************************************************************
* Generate a synthetic emission line (or ThAr) calibration spectrum as
* a series of equal height, equal width (in velocity space)
* Gaussians.
****************************************************************************/

#include <math.h>
#include "UVES_popler.h"
#include "gamm.h"
#include "const.h"
#include "error.h"

int UVES_combsynthThAr(cspectrum *cspec, params *par) {

  double   wlg=0.0,pref=0.0,sigl=0.0,sigr=0.0;
  int      fval=0,lval=0; /* First and last valid pixels in raw order */
  int      i=0;

  /* Find first valid pixel in raw array */
  i=0; while (cspec->st[i]<1) i++; fval=i;
  /* Find last valid pixel in raw array */
  i=cspec->np-1; while (cspec->st[i]<1) i--; lval=i;

  /* Loop over valid pixels */
  for (i=fval; i<=lval; i++) {
    /* Is this pixel near enough to a Gaussian for us to calculate the flux? */
    wlg=2.0*(double)((int)(cspec->wl[i]/2.0));
    wlg=(cspec->wl[i]-wlg>1.0) ? wlg+2.0 : wlg; 
    /* Initialise flux and error */
    cspec->fl[i]=cspec->no[i]=100.0; cspec->er[i]=cspec->ne[i]=0.1;
    cspec->co[i]=1.0;
      cspec->st[i]=1;
    if (C_C_K*fabs(cspec->wl[i]-wlg)/cspec->wl[i]<20.0) {
      pref=70000.0*C_FWHMSIG/wlg/C_SQRT2;
      sigl=(i) ? pref*(cspec->rwl[i-1]-wlg) : pref*(cspec->flwl-wlg);
      sigr=pref*(cspec->rwl[i]-wlg);
      cspec->fl[i]+=2000.0*(erffn(sigr)-erffn(sigl)); cspec->no[i]=cspec->fl[i];
      cspec->er[i]+=1.e-2*sqrt(cspec->fl[i]-100.0); cspec->ne[i]=cspec->er[i];
    }
    cspec->ef[i]=cspec->nf[i]=cspec->er[i];
  }

  return 1;

}
