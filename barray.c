/****************************************************************************
* Allocate an array of characters with length n. Function returns a
* pointer to first element
****************************************************************************/

#include <stdlib.h>
#include <stdbool.h>
#include "error.h"

bool *barray(unsigned long n) {
  
  bool   *ba=NULL;
  int      i=0;

  if (!(ba=(bool *)malloc((size_t)(n*sizeof(bool))))) {
    nferrormsg("barray(): Cannot allocate memory for array");
    return NULL;
  }

  /* Initialize array */
  for (i=0; i<n; i++) *(ba+i)=0.0;

  return ba;

}
