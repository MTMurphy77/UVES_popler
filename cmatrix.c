/****************************************************************************
* Allocate a matrix of characters with dimensions n x m. Function
* returns a pointer to (0,0) element
****************************************************************************/

#include <stdlib.h>
#include "memory.h"
#include "error.h"

char **cmatrix(unsigned long n, unsigned long m) {
  
  int      i=0;
  char     **cm=NULL;

  /* Allocate array of pointers */
  if (!(cm=(char **)malloc((size_t)(n*sizeof(char *))))) {
    nferrormsg("cmatrix(): Cannot allocate memory for pointer array");
    return NULL;
  }

  /* Allocate and initialize rows */
  if ((*cm=carray(n*m))==NULL) {
    nferrormsg("cmatrix(): Cannot allocate memory for array of size %dx%d",n,m);
    return NULL;
  }

  /* Set pointers to rows */
  for (i=1; i<n; i++) cm[i]=cm[i-1]+m;

  return cm;

}
