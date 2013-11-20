/***************************************************************************
SORT.H: Include file for sorting functions
***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */

/* STRUCTURES */
typedef struct TwoDble {
  double a;
  double b;
} twodble;

typedef struct DbleInt {
  double a;
  int    i;
} dbleint;

typedef struct DbleTwoInt {
  double a;
  int    i;
  int    j;
} dbletwoint;

/* FUNCTION PROTOTYPES */
double dselect(unsigned long k, unsigned long n, double *arr);
int qsort_absdarray(const void *x1, const void *x2);
int qsort_darray(const void *x1, const void *x2);
int qsort_dbleint(const void *dat1, const void *dat2);
int qsort_dbletwointarray(const void *dat1, const void *dat2);
int qsort_farray(const void *x1, const void *x2);
int qsort_twodarray(const void *dat1, const void *dat2);
int sort_2darray(int n, double *data1, double *data2);
