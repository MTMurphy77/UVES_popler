/****************************************************************************
* Qsort comparison routine for a single double precision array
****************************************************************************/

int qsort_darray(const void *x1, const void *x2) {

  if ((*(double *)x1) > (*(double *)x2)) return 1;
  else if ((*(double *)x1) == (*(double *)x2)) return 0;
  else return -1;

}
