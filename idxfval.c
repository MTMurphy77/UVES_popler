/****************************************************************************
* Find array index of first array element with value >= val;
****************************************************************************/

int idxfval(float *array, int n, float val) {

  int  i=-1;

  while (++i<n && val>array[i]);
  if (i<n) return i;
  return -1;

}
