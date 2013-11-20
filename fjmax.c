/****************************************************************************
* Find the maximum number in an array of floats
****************************************************************************/

#include "error.h"

float fjmax(float array[], int n) {

  int          i;
  float        max;

  if (n <= 0)
    errormsg("fjmax(): Bad arraysize passed.");

  max = array[0];
  for (i = 1; i < n; i++)
    if (array[i] > max)
      max = array[i];

  return max;

}
