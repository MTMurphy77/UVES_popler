/****************************************************************************
* Reject pixels from order edges which fall well below the general
* trend of errors determined from a median filtered version of the
* error array
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "memory.h"
#include "stats.h"
#include "error.h"

int UVES_order_rejsigedge(echorder *ord, params *par) {

  double   *med=NULL;
  int      nmed=0;
  int      i=0;

  /* Check first to see if order is useful */
  if (ord->nuse<MINUSE) return 1;

  /* Median filter length is specified as a fixed fraction of the order length */
  nmed=(MAX(NMEDERR,(1+2*(int)(par->ordmedfrac*0.5*(double)ord->np))));

  /* Allocate memory for median array */
  if ((med=darray(ord->np))==NULL)
    errormsg("UVES_order_rejsigedge(): Could not allocate memory\n\
\tfor med array of size %d",ord->np);

  /* Calculate median filtered error array */
  if (!medianrun(ord->er,med,ord->st,ord->np,nmed)) {
    nferrormsg("UVES_order_rejsigedge(): Error returned from medianrun()");
    return 0;
  }

  /* Reject pixels with errors well below median level */
  for (i=0; i<ord->np; i++)
    if (ord->st[i]==1 && ord->er[i]<med[i]*par->ordmedrej) ord->st[i]=OCLIP;

  /* Clean up */
  free(med);
  
  return 1;

}
