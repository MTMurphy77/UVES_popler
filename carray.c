/****************************************************************************
* Allocate an array of characters with length n. Function returns a
* pointer to first element
****************************************************************************/

#include <stdlib.h>
#include "error.h"

char *carray(unsigned long n) {
  
  int      i=0;
  char     *ca=NULL;

  if (!(ca=(char *)malloc((size_t)(n*sizeof(char))))) {
    nferrormsg("carray(): Cannot allocate memory for array");
    return NULL;
  }

  /* Initialize array */
  for (i=0; i<n; i++) *(ca+i)=' ';

  return ca;

}
