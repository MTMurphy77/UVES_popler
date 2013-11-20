/****************************************************************************
* Find the minimum number in an array of floats
****************************************************************************/

#include "error.h"

float fjmin(float array[], int n) {

  int          i;
  float        min;

  if (n <= 0)
    errormsg("fjmin(): Bad arraysize passed.");

  min = array[0];
  for (i = 1; i < n; i++)
    if (array[i] < min)
      min = array[i];

  return min;
}
