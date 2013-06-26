/****************************************************************************
* Find the initial continuum for all relevant orders
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "error.h"

int UVES_order_cont(spectrum *spec, int nspec, params *par) {

  int      i=0,j=0,k=0;

  /* Loop over spectra */
  for (i=0; i<nspec; i++) {
    /* Loop over useful orders */
    for (j=0; j<spec[i].nor; j++) {
      if (spec[i].or[j].nuse>=MINUSE) {
	/* Fit a continuum if order is useful */
	if (UVES_confit(spec[i].or[j].rdfl,spec[i].or[j].rder,
			spec[i].or[j].rdst,spec[i].or[j].nrdp,par->contftyp,
			par->contord,par->contsigl,par->contsigu,par->contpctl,
			1,spec[i].or[j].rdco)<1)
	  errormsg("UVES_order_cont(): Problem fitting cont. to order %d\n\
\tin file %s\n\t",j+1,spec[i].file);

	/* Go through order and recitfy pixels with no information */
	for (k=0; k<spec[i].or[j].nrdp; k++) {
	  if (spec[i].or[j].rdst[k]==RCLIP) {
	    spec[i].or[j].rdfl[k]=spec[i].or[j].rdco[k];
	    spec[i].or[j].rder[k]=-INFIN; spec[i].or[j].rdef[k]=-INFIN;
	  }
	}
      }
    }
  }
  
  return 1;
}
