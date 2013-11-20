/***************************************************************************
GAMM.H: Include file containing information for use with gamma, beta, and
error function calculations.
***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */
#define GAMM_ITMAX 1000    /* Maximum iterations for gammcf and gammser */
#define GAMM_EPS   3.0e-9  /* Relative precision */
#define GAMM_FPMIN 1.0e-30 /* Number near the smallest representable float */

/* STRUCTURES */

/* PROTOTYPES */
double betacf(double a, double b, double x);
double betai(double a, double b, double x);
double erfcc(double x);
double erffn(double x);
double erffc(double x);
double gammp(double a, double x);
double gammq(double a, double x);
double gammln(double x);
int gammcf(double *gcf, double a, double x, double *gln);
int gammser(double *gser, double a, double x, double *gln);
double gauss(double x);
double gauss_int(double A, double x0, double sigma, double a, double b);
