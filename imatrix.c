/****************************************************************************
* Allocate a matrix of integers with dimensions n x m. Function
* returns a pointer to (0,0) element
****************************************************************************/

#include <stdlib.h>
#include "memory.h"
#include "error.h"

int **imatrix(unsigned long n, unsigned long m) {
  
  int   **im=NULL;
  int   i=0;

  /* Allocate array of pointers */
  if (!(im=(int **)malloc((size_t)(n*sizeof(int *))))) {
    nferrormsg("imatrix(): Cannot allocate memory for pointer array");
    return NULL;
  }

  /* Allocate and initialize rows */
  if ((*im=iarray(n*m))==NULL) {
    nferrormsg("imatrix(): Cannot allocate memory for array of size %dx%d",n,m);
    return NULL;
  }

  /* Set pointers to rows */
  for (i=1; i<n; i++) im[i]=im[i-1]+m;

  return im;

}
