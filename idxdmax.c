/****************************************************************************
* Find the index of the maximum number in an array of n doubles. Returns -1
* if an error occurs, usually because an n<0 has been passed.
****************************************************************************/

int idxdmax(double *array, int n) {

  int          i=0,imax=0;
  double       max;

  if (n<0) return -1;
  max=array[0];
  for (i=1; i<n; i++) {
    if (array[i]>max) {
      max=array[i];
      imax=i;
    }
  }

  return imax;

}
