/****************************************************************************
* Qsort comparison routine for an array of double precision + integer
* structures
****************************************************************************/

#include "sort.h"

int qsort_dbleint(const void *dat1, const void *dat2) {

  if (((dbleint *)dat1)->a > ((dbleint *)dat2)->a) return 1;
  else if (((dbleint *)dat1)->a == ((dbleint *)dat2)->a) return 0;
  else return -1;

}
