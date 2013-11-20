/****************************************************************************
* Find array index of "first" array element with value >= val;
****************************************************************************/

int idxdval(double *array, int n, double val) {

  register int i=-1;

  while (++i<n && val>array[i]);
  if (i<n) return i;
  return -1;

}
