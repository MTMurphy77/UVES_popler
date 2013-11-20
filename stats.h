/***************************************************************************
STATS.H: Include file for the STATS function.
***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */
/* Functions */
#define MAX(a,b)   a>b ? a : b
#define MIN(a,b)   a<b ? a : b
#define MOD(a,b)   (b==0.0) ? b : a-((int)(a/b))*b

/* Parameters */
#define INFIN     1.e32 /* Effective infinity parameter                    */
#define KS_EPS1   1.e-3 /* Relative size of consecutuve series terms below
			   which ksprob() will stop calculating terms in
			   series                                          */
#define KS_EPS2   1.e-8 /* Relative total precision below which ksprob()
			   will stop calculating terms in series           */
#define KS_ITER   100   /* Maximum number of iterations allowed for
			   convergence in ksprob()                         */
#define MED_QR    0.341344746069 /* Value of quartile range for median()   */

/* STRUCTURES */
typedef struct Statset {
  double   mean;      /* Mean of sample                                    */
  double   rms;       /* RMS of sample                                     */
  double   prms;      /* Positive RMS of sample                            */
  double   nrms;      /* Negative RMS of sample                            */
  double   emean;     /* Error in mean of sample                           */
  double   wmean;     /* Weighted mean of sample                           */
  double   ewmean;    /* 1-sigma error in weighted mean                    */
  double   meansig;   /* Mean 1-sigma error in sample                      */
  double   eflwmean;  /* Expected fluctuation in weighted mean             */
  double   rmssig;    /* RMS of 1-sigma errors                             */
  double   chisq;     /* Chisquared about weighted mean                    */
  double   rchisq;    /* Reduced chisquared about weighted mean            */
  double   med;       /* Median of sample                                  */
  double   siqr;      /* 68% semi-interquartile range                      */
  double   emed;      /* Median of error array                             */
  double   esiqr;     /* 68% semi-interquartile range of errors array      */
  double   prun;      /* Probability from WW runs test                     */
} statset;

/* FUNCTION PROTOTYPES */
int autocorr_short(double *data, int n, double *corr, int len, int cntr);
int crank(unsigned long n, double *w, double *s);
int ftest(double *data1, unsigned long n1, double *data2, unsigned long n2,
	  double *f, double *prob);
double gasdev(long *idum);
int isodd(int number);
int kendall(double *data1, double *data2, unsigned long n, double *tau,
	    double *z, double *prob);
double ksprob(double alam);
int ksprob_2dist(double *data1, unsigned long n1, double *data2,
		 unsigned long n2, double *d, double *prob);
int ksprob_2fdist(float *data1, unsigned long n1, float *data2,
		  unsigned long n2, double *d, double *prob);
int median(double *dat, int *sts, int ndat, statset *stat, int opt);
int medianfilter(double *dat, double *err, int ndat, double medsig,
		 statset *stat, int *clip, int opt);
int medianrun(double *dat, double *med, int *sts, int ndat, int nfilt);
double ran(long *idum);
int ranarray(double *array, int *idx, int n, long *idum);
int sigclip(double *dat, double *err, double *efl, double *wgt, int ndat,
	    double clipsig, statset *stat, int *clip);
int spearman(double *data1, double *data2, unsigned long n, double *d, double *zd,
	     double *probd, double *rs, double *probrs);
int stats(double *dat, double *sig, double *efl, double *wgt, int *sts, int ndat,
	  int opt, statset *stat);
int utest(double *data1, unsigned long n1, double *data2, unsigned long n2,
	  double *U, double *prob);
int wwruns(double *dat, int *sts, int ndat, statset *stat, double med, int opt1,
	   int opt2);
