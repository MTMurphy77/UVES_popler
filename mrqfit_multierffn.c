/*************************************************************************** 
MRQFIT_MULTIERFFN: This is the fgauss algorithm taken from Numercial
Recipes but instead of just taking the value of a Gaussian at a given
x value, the Gaussian is intergrated between two x values using the
error function, then divided by the difference between the two x
values. This is usually to take into account the pixelization of a
Gaussian signal.  This particular algorithm fits mutliple gaussians
which are all related to a primary (e.g. doublets). Note using a
completely unconstrained fit is no different from fitting each
absorption in the multiplet separately (as in mrqfit_erffn.c).
Parameters of daughter Gaussians are ratio's to the primary
gaussian. Constraints should be placed as such and are passed to the
algorithm by lim[][2]. At the moment you can not place limits on
constraints (i.e. a doublet ratio limited between 1 and 2) - This
could probably be implemented by use of an additional two na x na
matrices added onto lim[][].  Here's what NR has to say about it:

"y(x,a)is the sum of na/3 Gaussians (15.5.16). The amplitude, center, and
width of the Gaussians are stored in consecutive locations of a: a[i] = Bk,
a[i+1] = Ek, a[i+2] = Gk, k = 1, ..., na/3. The dimensions of the arrays
are a[1..na], dyda[1..na]."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.
2) Have made the first argument an array. However, for a gaussian fit,
   only the first element should be filled and nx should be set to 1.
3) Have added two more parameters per gaussian: 4th parameter is the
   continuum slope and the 5th parameter is the continuum
   level. Therefore, the user must make sure there are na/5 gaussians
   being fit.
4) To clarify, the parameters in a are defined with respect to the
   following equation: 
   y(x) = a0*exp[-(x-a1)^2/2*a2^2] + a3*(x-a1) + a4 +
          a0*a5*exp[-((x-a1)*a6)^2/2*(a2*a7)^2] + a3*a8*(x-a1)*a6 + a4*a9 + ...
5) Have added an na*3 matrix additional argument defining the limits of the
   parameters and constraints on daughter Gaussians.
6) If the dyda array is passed as NULL then no attempt is made to fill
   it with the derivatives.
7) Spheres of influence of sets of parameters can be specified for
   each set of parameters (i.e. each Gaussian). If the spheres of
   influence overlap then for the Gaussians the sum is taken, but for
   the continua an average is taken. Sphere of influence for each
   parameter set should be specified in the matrix entry corresponding
   to the first parameter of that set. If you wish to always include a
   particular parameter set then either set its bounds to -1. If a
   NULL **inf matrix is passed it is assumed that all sets contribute
   to all fits. All pixels are zero-indexed.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "fit.h"
#include "gamm.h"
#include "memory.h"
#include "const.h"
#include "error.h"

int mrqfit_multierffn(double *x, double *a, double **lim, int *ia, int **inf,
		      double *y, double *dyda, int nx, int na, int idx) {

  double bma=0.0,bsmas=0.0,aexp=0.0,bexp=0.0,aexps=0.0,bexps=0.0;
  double expa=0.0,expb=0.0,erfa=0.0,erfb=0.0;
  double *adum=NULL;
  int    *inc=NULL,ninc=0,nset=0;
  int    i=0,j=0,k=0;

  /* Initialise adum */
  if ((adum=darray(na))==NULL)
    errormsg("mrqfit_multierffn(): Cannot allocate memory for adum\n\
\tarray of size %d",na);

  /* Check to make sure 2 x limits have been entered */
  if (nx!=2) {
    nferrormsg("mrqfit_multierffn(): Only %d x-limits provided.\n\
\tMust provide 2.",nx); return 0;
  }

  /* Check to make sure number of parameters is divisible by 5 */
  if ((MOD(na,5))) {
    nferrormsg("mrqfit_multierffn(): Number of input parameters (=%d)\n\
\tnot a multiple of 5",na); return 0;
  }

  /* Check to make sure x-limits are not the same */
  bma=x[1]-x[0]; bsmas=x[1]*x[1]-x[0]*x[0];
  if (x[0]==x[1]) {
    nferrormsg("mrqfit_multierffn(): x-limits are equal (=%lf)",x[0]);
    return 0;
  }

  /* Check to make sure primary is not constrained with respect to itself */
  for (i=0;i<5;i++) {
    if (lim[i][2] && na>5) 
      errormsg("mrqfit_multierffn(): Cannot constrain primary gaussian with\n\
respect to itself. Set lim[%d][2] = 0",i);
  }

  /* Ensure if a daughter Gaussian is constrained that it is not
     allowed to vary */
  for (i=5;i<na;i+=5) {
    for (j=0;j<5;j++) {
      if (lim[i+j][2] && ia[i+j])
	errormsg("mrqfit_multierffn(): Nonsensical restraints placed on\n\
\tparameter a[%d]",i+j);
    }
  }

  /** Make some checks on the spheres of influence of each parameter set **/
  /* Initialise inc */
  ninc=0;
  nset=na/5;
  if ((inc=iarray(nset))==NULL)
    errormsg("mrqfit_multierffn(): Cannot allocate memory for inc\n\
\tarray of size %d",nset);
  /* Set inc */
  if (inf!=NULL) {
    for (i=0; i<nset; i++) {
      if ((idx>=inf[5*i][0] && idx<=inf[5*i][1]) || inf[5*i][0]==-1) {
	inc[i]=1; ninc++;
      }
      else inc[i]=0;
    }
  }
  else for (i=0; i<nset; i++) { inc[i]=1; ninc++; }

  /** Add up contributions from all Gaussians **/
 
  /* Force parameters to lie within limits */
  for (j=0; j<na; j++) {
    if (ia[j]) { a[j]=(MAX(a[j],lim[j][0])); a[j]=(MIN(a[j],lim[j][1])); }
  }

  /* Primary Gaussian independantly, first */
  i=0; *y=0.0;
  /* Check to make sure width of Gaussian isn't zero */
  if (fabs(a[2])<FPMIN) {
    nferrormsg("mrqfit_multierffn(): Gaussian %d has near zero\n\
\tor negative width (=%lf)",1,a[2]); return 0;
  }

  /* Calculate function and it's derivatives */
  aexp=(x[0]-a[1])/(C_SQRT2*a[2]); bexp=(x[1]-a[1])/(C_SQRT2*a[2]);
  aexps=aexp*aexp; bexps=bexp*bexp;
  expa=(aexps<80.0) ? exp(-aexps) : 0.0; expb=(bexps<80.0) ? exp(-bexps) : 0.0;
  erfa=C_SQRTPI/C_SQRT2*erffn(aexp); erfb=C_SQRTPI/C_SQRT2*erffn(bexp);
  *y += (double)inc[0]*a[0]*a[2]*(erfb-erfa)/bma +
    (double)inc[0]/(double)ninc*(a[3]/2.0*bsmas/bma - a[3]*a[1] + a[4]);

  /* Note: derivatives incomplete and depend on daughters */
  if (dyda!=NULL) {
    dyda[0]=a[2]*(erfb-erfa)/bma;
    dyda[1]=a[0]*(expa-expb)/bma - a[3];
    dyda[2]=a[0]/bma*(C_SQRT2*(aexp*expa-bexp*expb) + (erfb-erfa));
    dyda[3]=bsmas/2.0/bma - a[1];
    dyda[4]=1.0;
  }

  /* Daughter Gaussians */
  for (i=5,k=1; i<na-1; i+=5,k++) {
    /* Check to make sure width of Gaussian isn't zero */
    if (fabs(a[2]*a[i+2])<FPMIN) {
	nferrormsg("mrqfit_multierffn(): Gaussian %d has near zero\n\
\tor negative width (=%lf)",i/5+1,a[2]*a[i+2]); return 0;
    }

    /* Create adum to deal implementing a variable number of constraints */
    for (j=0;j<5;j++) {
      if (lim[i+j][2]) adum[i+j] = a[j]*a[i+j];
      else adum[i+j] = a[i+j];
    }

    /* Calculate function and it's derivatives */
    /* If constrained then there is no derivative for that param, but
       there is an affect on the primary derivative due to a[0+j].  If
       unconstrained then no effect on primary derivatives exist. */

    aexp=(x[0]-adum[i+1])/(C_SQRT2*adum[i+2]);
    bexp=(x[1]-adum[i+1])/(C_SQRT2*adum[i+2]);

    aexps=aexp*aexp; bexps=bexp*bexp;
    expa=(aexps<80.0) ? exp(-aexps) : 0.0; expb=(bexps<80.0) ? exp(-bexps) : 0.0;
    erfa=C_SQRTPI/C_SQRT2*erffn(aexp); erfb=C_SQRTPI/C_SQRT2*erffn(bexp);

    *y += (double)inc[k]*adum[i]*adum[i+2]*(erfb-erfa)/bma + 
      (double)inc[k]/(double)ninc*
      (adum[i+3]/2.0*bsmas/bma - adum[i+3]*adum[i+1] + adum[i+4]);

    if (dyda!=NULL) {
      if (lim[i][2] && ia[0]) dyda[0]+=a[i]*adum[i+2]*(erfb-erfa)/bma;
      if (lim[i+1][2] && ia[1]) dyda[1]+=adum[i]*(expa-expb)/bma - adum[i+3]*a[i+1];
      if (lim[i+2][2] && ia[2]) dyda[2]+=adum[i]*a[i+2]/bma*
	((erfb-erfa) + C_SQRT2*(aexp*expa-bexp*expb));
      if (lim[i+3][2] && ia[3]) dyda[3]+=a[i+3]/2.0*bsmas/bma - a[i+3]*adum[i+1];
      if (lim[i+4][2] && ia[4]) dyda[4]+=a[i+4];
      dyda[i] = (!lim[i][2]) ? adum[i+2]*(erfb-erfa)/bma : 0.0;
      dyda[i+1] = (!lim[i+1][2]) ? adum[i]*(expa-expb)/bma - adum[i+3] : 0.0;
      dyda[i+2] = (!lim[i+2][2]) ? adum[i]/bma*
	((erfb-erfa) + C_SQRT2*(aexp*expa - bexp*expb)) : 0.0;
      dyda[i+3] = (!lim[i+3][2]) ? bsmas/bma/2.0 - adum[i+1] : 0.0;
      dyda[i+4] = (!lim[i+4][2]) ? 1.0 : 0.0;
    }
  }
  free(adum); free(inc);
  return 1; 
}
