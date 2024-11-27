/***************************************************************************
MEMORY.H: Include file for allocation of memory routines.
***************************************************************************/

/* INCLUDE FILES */
#include <stdbool.h>

/* DEFINITIONS */

/* STRUCTURES */
typedef struct DblePair{
  double x1;
  double x2;
} dblepair;

/* PROTOTYPES */
bool *barray(unsigned long n);
char *carray(unsigned long n);
char **cmatrix(unsigned long n, unsigned long m);
double *darray(unsigned long n);
double **dmatrix(unsigned long n, unsigned long m);
double ***dcube(unsigned long n, unsigned long m, unsigned long p);
double ****d4cube(unsigned long n, unsigned long m, unsigned long p, unsigned long q);
float *farray(unsigned long n);
float **fmatrix(unsigned long n, unsigned long m);
float ***fcube(unsigned long n, unsigned long m, unsigned long p);
int *iarray(unsigned long n);
int **imatrix(unsigned long n, unsigned long m);
int ***icube(unsigned long n, unsigned long m, unsigned long p);
