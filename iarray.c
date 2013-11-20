/****************************************************************************
* Allocate an array of integers with length n. Function returns a pointer
* to first element
****************************************************************************/

#include <stdlib.h>
#include "error.h"

int *iarray(unsigned long n) {
  
  int   *ia=NULL;
  int   i=0;

  if (!(ia=(int *)malloc((size_t)(n*sizeof(int))))) {
    nferrormsg("iarray(): Cannot allocate memory for array");
    return NULL;
  }

  /* Initialize array */
  for (i=0; i<n; i++) *(ia+i)=0;

  return ia;

}
