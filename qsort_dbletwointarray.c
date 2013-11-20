/****************************************************************************
* Qsort comparison routine for an array of double precision + 2 integer
* structures
****************************************************************************/

#include "sort.h"

int qsort_dbletwointarray(const void *dat1, const void *dat2) {

  if (((dbletwoint *)dat1)->a > ((dbletwoint *)dat2)->a) return 1;
  else if (((dbletwoint *)dat1)->a == ((dbletwoint *)dat2)->a) return 0;
  else return -1;

}
