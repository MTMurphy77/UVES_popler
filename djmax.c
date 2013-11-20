/****************************************************************************
* Find the maximum number in an array of doubles
****************************************************************************/

#include "error.h"

double djmax(double array[], int n) {

  int          i;
  double       max;

  if (n <= 0)
    errormsg("djmax(): Bad arraysize passed.");

  max = array[0];
  for (i = 1; i < n; i++)
    if (array[i] > max)
      max = array[i];

  return max;
}
