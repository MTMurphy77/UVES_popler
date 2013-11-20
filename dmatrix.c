/****************************************************************************
* Allocate a matrix of double precision with dimensions n x m. Function
* returns a pointer to (0,0) element
****************************************************************************/

#include <stdlib.h>
#include "memory.h"
#include "error.h"

double **dmatrix(unsigned long n, unsigned long m) {
  
  double   **dm=NULL;
  int      i=0;

  /* Allocate array of pointers */
  if (!(dm=(double **)malloc((size_t)(n*sizeof(double *))))) {
    nferrormsg("dmatrix(): Cannot allocate memory for pointer array");
    return NULL;
  }

  /* Allocate and initialize rows */
  if ((*dm=darray(n*m))==NULL) {
    nferrormsg("dmatrix(): Cannot allocate memory for array of size %dx%d",n,m);
    return NULL;
  }

  /* Set pointers to rows */
  for (i=1; i<n; i++) dm[i]=dm[i-1]+m;

  return dm;

}
