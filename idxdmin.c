/****************************************************************************
* Find the index of the minimum number in an array of n doubles. Returns -1
* if an error occurs, usually because an n<0 has been passed.
****************************************************************************/

int idxdmin(double *array, int n) {

  int          i=0,imin=0;
  double       min;

  if (n<0) return -1;
  min=array[0];
  for (i=1; i<n; i++) {
    if (array[i]<min) {
      min=array[i];
      imin=i;
    }
  }

  return imin;

}
