/****************************************************************************
* Apply scalings of flux and error arrays, either before (opt=0) or
* after (opt=1) orders/spectra are combined
****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "UVES_popler.h"
#include "error.h"

int UVES_scale(spectrum *spec, int opt) {

  double scale=0.0;
  int    i=0,j=0;
  int    flxerr=0,befaft=0;
  char   insclstr[VHUGESTRLEN]="\0";
  char   *token;

  /* Initialize order structures for storing scaling information and
     fill those structures from the scaling string in the spectrum
     structure */
  if (!opt) {
    /* If this is pre-combination, first task is to initialise the
       scaling information for each order */
    for (i=0; i<spec->nor; i++) {
      spec->or[i].insclfe=spec->or[i].insclbefaft=0;
      spec->or[i].inscl=1.0;
    }
    /* Read scaling string in the spectrum structure and fill
       information in order structure for those orders where scaling
       is required */
    if (spec->inscl[0]!='\0') {
      sprintf(insclstr,"%s",spec->inscl);
      j=0; token=strtok(insclstr,";");
      while (token!=NULL) {
	j++;
	if (sscanf(token,"%d,%d,%lf,%d",&i,&flxerr,&scale,&befaft)!=4) {
	  nferrormsg("UVES_scale(): Error reading section %d of flux/error scaling\n\
\tstring, '%s'. Format needs to be '[int],[int],[dbl],[int]'.",j,token);
	  return 0;
	}
	if (i<-1 || i==0 || i>spec->nor) {
	  nferrormsg("UVES_scale(): Echelle order specified, %d, in section\n\
\t%d of flux/error scaling string, '%s', is out of range.\n\
\tSpectrum has %d orders.",i,j,token); return 0;
	}
	if (i!=-1) {
	  i--; spec->or[i].insclfe=flxerr; spec->or[i].insclbefaft=befaft;
	  spec->or[i].inscl=scale;
	} else {
	  for (i=0; i<spec->nor; i++) {
	    spec->or[i].insclfe=flxerr; spec->or[i].insclbefaft=befaft;
	    spec->or[i].inscl=scale;
	  }
	}
	token=strtok(NULL,";");
      }
    }
  }

  /* Now do the scaling itself, scaling the original or redispersed
     arrays depending on whether this is before or after
     combination */
  if (!opt) {
    /* Scale the original arrays */
    for (i=0; i<spec->nor; i++) {
      if (spec->or[i].insclbefaft==1) {
	if (spec->or[i].insclfe==0 || spec->or[i].insclfe==2) {
	  for (j=0; j<spec->or[i].np; j++)
	    spec->or[i].fl[j]*=spec->or[i].inscl;
	}
	if (spec->or[i].insclfe==1 || spec->or[i].insclfe==2) {
	  for (j=0; j<spec->or[i].np; j++)
	    spec->or[i].er[j]*=spec->or[i].inscl;
	}
      }
    }
  } else {
    /* For post-combination scaling */
    for (i=0; i<spec->nor; i++) {
      if (spec->or[i].insclbefaft==2) {
	if (spec->or[i].insclfe==0 || spec->or[i].insclfe==2) {
	  for (j=0; j<spec->or[i].nrdp; j++)
	    spec->or[i].rdfl[j]*=spec->or[i].inscl;
	}
	if (spec->or[i].insclfe==1 || spec->or[i].insclfe==2) {
	  for (j=0; j<spec->or[i].nrdp; j++) {
	    spec->or[i].rder[j]*=spec->or[i].inscl;
	    spec->or[i].rdef[j]*=spec->or[i].inscl;
	    spec->or[i].rdme[j]*=spec->or[i].inscl;
	  }
	}
      }
    }
  }

  return 1;

}
