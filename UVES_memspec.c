/****************************************************************************
* Allocate or free memory for raw arrays for all echelle orders in a specific
* spectrum
*
* Opt=i allocates memory to echelle order with index i (i.e. order number i+1).
* Opt=[Number of echelle orders in spectrum] allocates memory and initializes
*     status arrays to 1 for all echelle orders in the spectrum.
* Opt=-i frees memory from echelle order number i (i.e. with index i-1).
* Opt=-1-[Number of echelle orders] frees memory from all echelle orders
*     in the spectrum.
*
* NOTE: When freeing memory, this routine only checks whether MINUSE
* valid pixels are present in a given echelle order when freeing
* memory from all echelle orders of a spectrum. When freeing memory
* from a single, user-specified order, it is assumed that the user has
* already checked that the memory is there to be freed before calling
* this routine.
****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

#define ERR_ARRAY { \
    nferrormsg("UVES_memspec(): Cannot allocate memory for %s\n\
\tarray of order %d of size %d for file\n\t%s",array_name,i+1,spec->or[i].np,spec->file); \
    return 0; }

int UVES_memspec(spectrum *spec, params *par, int opt) {

  int      os=0,oe=0;
  int      i=0,j=0;
  char     array_name[NAMELEN]="\0";

  /* If memory is to be allocated */
  if (opt>=0 && opt<=spec->nor) {
    os=(opt<spec->nor) ? opt : 0;
    oe=(opt<spec->nor) ? opt : spec->nor-1;
    /* Loop over nominated echelle orders */
    for (i=os; i<=oe; i++) {
      /* Allocate memory for raw data arrays for this echelle order */
      /* Unlike most other instruments/pipelines, HIRES REDUX and
	 MAGE files have an actual wavelength array to be read in,
	 so allocate memory for the raw wavelength array first */
      if (spec->ftype==FTHIRX || spec->ftype==FTMAGE) {
	sprintf(array_name,"%s","raw wavel.");
	if ((spec->or[i].wl=darray(spec->or[i].np))==NULL) { ERR_ARRAY; }
      }    
      sprintf(array_name,"%s","vac. wavel.");
      if ((spec->or[i].vhwl=darray(spec->or[i].np))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"%s","right-edge vac. wavel.");
      if ((spec->or[i].vhrwl=darray(spec->or[i].np))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"%s","flux");
      if ((spec->or[i].fl=darray(spec->or[i].np))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"%s","error");
      if ((spec->or[i].er=darray(spec->or[i].np))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"%s","resolution");
      if ((spec->or[i].res=darray(spec->or[i].np))==NULL) { ERR_ARRAY; }
      if (par->thar==1) {
	sprintf(array_name,"%s","ThAr flux");
	if ((spec->or[i].th=darray(spec->or[i].np))==NULL) { ERR_ARRAY; }
	sprintf(array_name,"%s","ThAr error");
	if ((spec->or[i].ter=darray(spec->or[i].np))==NULL) { ERR_ARRAY; }
      }
      sprintf(array_name,"%s","status");
      if ((spec->or[i].st=iarray(spec->or[i].np))==NULL) { ERR_ARRAY; }
      if (par->thar==1) {
	sprintf(array_name,"%s","ThAr status");
	if ((spec->or[i].tst=iarray(spec->or[i].np))==NULL) { ERR_ARRAY; }
      }
      /* Initialise status array(s) */
      for (j=0; j<spec->or[i].np; j++) spec->or[i].st[j]=1;
      if (par->thar==1) for (j=0; j<spec->or[i].np; j++) spec->or[i].tst[j]=1;
    }
  } else if (opt<0 && opt>=-1-spec->nor) {
    os=(opt>-1-spec->nor) ? -1-opt : 0;
    oe=(opt>-1-spec->nor) ? -1-opt : spec->nor-1;
    /* Loop over nominated echelle orders */
    for (i=os; i<=oe; i++) {
      if (os==oe || spec->or[i].nuse>=MINUSE) {
	/* If memory is to be freed */
	if (spec->ftype==FTHIRX || spec->ftype==FTMAGE) free(spec->or[i].wl);
	free(spec->or[i].vhwl); free(spec->or[i].vhrwl); free(spec->or[i].fl);
	free(spec->or[i].er); free(spec->or[i].res); free(spec->or[i].st);
	if (par->thar==1) {
	  free(spec->or[i].th); free(spec->or[i].ter); free(spec->or[i].tst);
	}
      }
    }
  }

  return 1;

}
