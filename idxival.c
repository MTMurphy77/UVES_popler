/****************************************************************************
* Find array index of val; returns -1 if val not found
****************************************************************************/

int idxival(int array[], int n, int val) {

  int  i;

  for (i = 0; i < n; i++)
    if (array[i] == val)
      return i;

  return -1;
}
