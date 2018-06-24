/****************************************************************************
* Do some simple statistics on a spectrum order
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "stats.h"
#include "memory.h"
#include "error.h"

int UVES_order_stats(echorder *ord, params *par) {

  double   *data=NULL;    /* Array on which statistics are determined */
  int      fval=0,lval=0,nval=0;
  int      i=0;
  statset  stat;

  /* Check first to see if order is useful */
  if (ord->nuse<MINUSE) return 1;

  /* Allocate memory to data array */
  if ((data=darray(ord->nrdp))==NULL)
    errormsg("UVES_order_stats(): Could not allocate memory\n\
\tfor data array of size %d",ord->nrdp);

  /* Find first and last valid pixels in order */
  i=0; while (ord->rdst[i]<1) i++; fval=i;
  i=ord->nrdp-1; while (ord->rdst[i]<1) i--; lval=i;

  /* Load the data array with valid pixels */
  for (nval=0,i=fval; i<=lval; i++) {
    if (ord->rdst[i]==1 && ord->rder[i]>DRNDTOL) {
      data[nval++]=ord->rdfl[i]/ord->rder[i];
    }
  }
  /* Find median S/N of valid pixels in order */
  if (!median(data,NULL,nval,&stat,0)) {
    nferrormsg("UVES_order_stats(): Error returned from median()"); return 0;
  }
  ord->medsnr=stat.med;

  /* Calculate a median-smoothed error array */
  if (!medianrun(ord->rder,ord->rdme,ord->rdst,ord->nrdp,NMEDERR)) {
    nferrormsg("UVES_order_stats(): Error returned from medianrun()");
    return 0;
  }
  /* Calculate the median-smoothed alternative error array for
     pre-v0.74 backwards compatibility */
  if (par->version<0.74) {
    if (!medianrun(ord->rder_a,ord->rdme_a,ord->rdst,ord->nrdp,NMEDERR)) {
      nferrormsg("UVES_order_stats(): Error returned from medianrun()");
      return 0;
    }
  }

  /* Clean up */
  free(data);
  
  return 1;

}
