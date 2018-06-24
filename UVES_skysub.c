/****************************************************************************
* Sky-subtract all orders in a spectrum.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

int UVES_skysub(spectrum *spec, cspectrum *cspec, params *par) {

  double   tol=2.e-4,swll=0.0,swlm=0.0,ewll=0.0,ewlm=0.0;
  int      i=0,j=0,k=0,l=0,m=0;

  /* Check first that the spectrum actually needs sky-subtracting */
  if (spec->skysub) {
    /* Loop over all orders in OBJECT spectrum */
    for (l=0; l<spec->nor; l++) {
      /* Check first to see if object order is useful */
      if (spec->or[l].nuse>=MINUSE) {
	/** It is NOT presumed that each order in the sky spectrum
	    really matches the corresponding one in the object
	    spectrum. Therefore, we need to find a matching sky
	    echelle order based on its wavelength coverage compared to
	    each object order. This makes use of a hard-coded relative
	    tolerance level, tol, for both the starting and ending
	    wavelengths of the order **/
	swll=cspec->wl[spec->or[l].csidx]; ewll=cspec->wl[spec->or[l].ceidx];
	/* Find a matching sky order */
	for (m=0; m<spec->ss.nor; m++) {
	  /* Make sure that sky order is useful */
	  if (spec->ss.or[m].nuse>=MINUSE) {
	    swlm=cspec->wl[spec->ss.or[m].csidx]; ewlm=cspec->wl[spec->ss.or[m].ceidx];
	    /* Check if starting and ending order wavelengths are
	       within tolerance of object order's */
	    if (fabs((swlm-swll)/swll)<tol && fabs((ewlm-ewll)/ewll)<tol) break;
	  }
	}
	/* If a match was found, continue with sky subtraction */
	if (m<spec->ss.nor) {
	  /* Loop over redispersed pixels in OBJECT spectrum */
	  for (i=0,j=spec->or[l].csidx; i<spec->or[l].nrdp; i++,j++) {
	    k=j-spec->ss.or[m].csidx;
	    /* If there's no sky pixel corresponding to this object
	       pixel then set the status array for the object pixel
	       appropriately */
	    if (j<spec->ss.or[m].csidx || j>spec->ss.or[m].ceidx)
	      spec->or[l].rdst[i]=RCLIP;
	    /* Otherwise, check to make sure both the sky pixel status
	       is valid; if not, set the status of the object pixel
	       appropriately */
	    else if (spec->ss.or[m].rdst[k]!=1) spec->or[l].rdst[i]=RCLIP;
	    /* Otherwise, if the object pixel is valid, do the sky
	       subtraction */
	    else if (spec->or[l].rdst[i]==1) {
	      spec->or[l].rdfl[i]-=spec->ss.or[m].rdfl[k];
	      spec->or[l].rder[i]=sqrt(spec->or[l].rder[i]*spec->or[l].rder[i]+
				       spec->ss.or[m].rder[k]*spec->ss.or[m].rder[k]);
	      spec->or[l].rdef[i]=sqrt(spec->or[l].rdef[i]*spec->or[l].rdef[i]+
				       spec->ss.or[m].rdef[k]*spec->ss.or[m].rdef[k]);
	      /* Do the subtraction for alternative arrays for
		 pre-v0.74 backwards compatibility */
	      spec->or[l].rdfl_a[i]-=spec->ss.or[m].rdfl[k];
	      spec->or[l].rder_a[i]=sqrt(spec->or[l].rder_a[i]*spec->or[l].rder_a[i]+
					 spec->ss.or[m].rder[k]*spec->ss.or[m].rder[k]);
	      spec->or[l].rdef_a[i]=sqrt(spec->or[l].rdef_a[i]*spec->or[l].rdef_a[i]+
					 spec->ss.or[m].rdef[k]*spec->ss.or[m].rdef[k]);
	    }
	  }
	} else {
	  /* If no matching sky order was found, set the number of
	     useful pixels in it to zero and set the status array of
	     all pixels in obejct order appropriately */
	  spec->or[l].nuse=0;
	  for (i=0; i<spec->or[l].nrdp; i++) spec->or[l].rdst[i]=RCLIP;
	}
      }
    }
  }

  return 1;

}
