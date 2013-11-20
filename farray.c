/****************************************************************************
* Allocate an array of double precision with length n. Function returns a
* pointer to first element
****************************************************************************/

#include <stdlib.h>
#include "error.h"

float *farray(unsigned long n) {
  
  float    *fa=NULL;
  int      i=0;

  if (!(fa=(float *)malloc((size_t)(n*sizeof(float))))) {
    nferrormsg("farray(): Cannot allocate memory for array");
    return NULL;
  }

  /* Initialize array */
  for (i=0; i<n; i++) *(fa+i)=0.0;

  return fa;

}
