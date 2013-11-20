/****************************************************************************
* Find the maximum number in an array of integers
****************************************************************************/

#include "error.h"

int ijmax(int array[], int n) {

  int          i;
  int          max;

  if (n <= 0)
    errormsg("ijmax(): Bad arraysize: %d", n);

  max = array[0];
  for (i = 1; i < n; i++)
    if (array[i] > max)
      max = array[i];

  return max;
}
