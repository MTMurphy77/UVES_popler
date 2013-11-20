/****************************************************************************
* Allocate an array of double precision with length n. Function returns a
* pointer to first element
****************************************************************************/

#include <stdlib.h>
#include "error.h"

double *darray(unsigned long n) {
  
  double   *da=NULL;
  int      i=0;

  if (!(da=(double *)malloc((size_t)(n*sizeof(double))))) {
    nferrormsg("darray(): Cannot allocate memory for array");
    return NULL;
  }

  /* Initialize array */
  for (i=0; i<n; i++) *(da+i)=0.0;

  return da;

}
