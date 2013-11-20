/****************************************************************************
* Find the minimum number in an array of doubles
****************************************************************************/

#include "error.h"

double djmin(double array[], int n) {

  int          i;
  double       min;

  if (n <= 0)
    errormsg("djmin(): Bad arraysize passed.");

  min = array[0];
  for (i = 1; i < n; i++)
    if (array[i] < min)
      min = array[i];

  return min;
}
