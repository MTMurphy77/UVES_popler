/***************************************************************************
FIT.H: Include file for various fitting functions and routines
***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */
/* Functions */
#define MAX(a,b)   a>b ? a : b
#define MIN(a,b)   a<b ? a : b
#define MOD(a,b) (b==0.0) ? b : a-((int)(a/b))*b
#define SQR(a) ((a) == 0.0 ? 0.0 : (a)*(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)

/* Parameters */
#define MRQNMAX     100       /* Maximum number of mrqmin iterations */ 
#define SVDNMAX     200       /* Maximum number of svdfit iterations */
#define BRNMAX      100       /* Maximum number of brent iterations  */
#define DBNMAX      100       /* Maximum number of dbrent iterations */
#define ZBNMAX      100       /* Maximum number of zbrent iterations */
#define NBTOL         3       /* Bad pixel threshold for EW.c doublet fits */
#define FPMIN         1.e-30  /* Number near smallest representable float */ 
#define BR_ACC        1.e-12  /* brent fractional accuracy */
#define DB_ACC        1.e-12  /* dbrent fractional accuracy */
#define ZB_EPS        1.e-30  /* zbrent machine precision */
#define FXY_ACC       1.e-10  /* fitexy accuracy */
#define SVDFIT_TOL    1.e-20  /* Tolerance level for svdfit() */
#define FIT_INFIN     1.e32   /* Effective infinity */
#define MNB_GLIMIT  100.0     /* Maximum magnification allowed for a parabolic-fit */

/* Flags */
#define GJSNGLR      -1       /* Flag denoting that a singular gaussj() matrix 
				 was formed */
#define MRQITMAX     -2       /* Flag indicating that over MRQNMAX
				 iterations were required for fit */

/* STRUCTURES */
typedef struct FitSet {
  double *x;           /* x data array                                      */
  double *sx;          /* x sigma data array                                */
  double *wx;          /* x weights array                                   */
  double *y;           /* y data array                                      */
  double *sy;          /* y sigma data array                                */
  double *wy;          /* y weights array                                   */
  double *p;           /* Polynomial fit coefficients (or other useful no.s)*/
  int    n;            /* Number of elements in data arrays                 */
  int    ns;           /* Number of valid (s=1) elements in data arrays     */
  int    nc;           /* Number of unclipped elements in clipped arrays    */
  int    *s;           /* Status array (s=1 for valid elements)             */
  int    *c;           /* Clip array (c=1 for unclipped elements)           */
} fitset;

/* FUNCTION PROTOTYPES */
int brent(fitset *fit, double ax, double bx, double cx, double (*f)(fitset *,double),
	  double tol, double *xmin, double *soln);
double chebyshev_eval(double x, double *coef, int ord);
double chixy (fitset *fit, double bang);
double chixy_fixa (fitset *fit, double bang);
int covsrt(double **covar, int ma, int *ia, int mfit);
int dpolcof(double *xa, double *ya, int n, double *cof);
int dpolint(double *xa, double *ya, int n, double x, double *y, double *dy);
int fitexy(double *x, double *y, double *sigx, double *sigy, int *stsx, int *stsy,
	   int ndat, double sigclip, int nclip, double *a, double *b, double *siga,
	   double *sigb, double *chi2, int opt1, int opt2, int verb);
int gaussj(double **a, int n, double **b, int m);
double legendre_eval(double x, double *coef, int ord);
int linfit(double *x, double *y, double *sig, int *sts, int ndata, int opt1, int opt2,
	   double *a, double *b, double *siga, double *sigb, double *chi2);
int linreg(double *x, double *y, int n, double *a, double *b);
int mnbrak(fitset *fit, double *ax, double *bx, double *cx, double *fa, double *fb,
	   double *fc, double (*func)(fitset *, double));
int mrqfit_erffn(double *x, double *a, double **lim, int *ia, int **inf,
		 double *y, double *dyda, int nx, int na, int idx);
int mrqfit_gauss(double *x, double *a, double **lim, int *ia, int **inf,
		 double *y, double *dyda, int nx, int na, int idx);
int mrqfit_multierffn(double *x, double *a, double **lim, int *ia, int **inf,
		      double *y, double *dyda, int nx, int na, int idx);
int mrqcof(double **x, double *y, double *sig, int ndata, int nx, double *a,
	   double **lim, int *ia, int **inf, int ma, double **alpha,
	   double *beta, double *chisq,
	   int (*funcs)(double *, double *, double **, int *, int **, double *,
			double *, int, int, int));
int mrqmin(double **x, double *y, double *sig, int ndata, int nx, double *a,
	   double **lim, int *ia, int **inf, int ma, double **covar, 
	   double **alpha, double *chisq, double *alamda,
	   int (*funcs)(double *, double *, double **, int *, int **, double *,
			double *, int, int, int));
double poly_eval(double x, double *coef, int ord);
double poly_nooffset_eval(double x, double *coef, int ord);
double pythag(double a, double b);
int spline(double *x, double *y, int n, double yp1, double ypn, double *y2,
	   int opt);
double splint(double *xa, double *ya, double *y2a, int n, double x);
int splint_array(double *xa, double *ya, double *y2a, int na, double *x, double *y,
		 int n);
int svbksb(double **u, double *w, double **v, int m, int n, double *b,
	   double *x);
int svdcmp(double **a, int m, int n, double *w, double **v);
int svdfit(double *x, double *y, double *sig, int ndata, double *a, int ma,
	double ***u, double ***v, double **w, double *chisq,
	int (*funcs)(double, double *, int));
int svdfit_chebyshev(double x, double *pt, int nt);
int svdfit_legendre(double x, double *pl, int nl);
int svdfit_poly(double x, double *p, int np);
int svdfit_poly_nooffset(double x, double *p, int np);
int svdvar(double **v, int ma, double *w, double ***cvm);
int zbrent(fitset *fit, double (*func)(fitset *, double), double x1, double x2,
	   double tol, double *root);
