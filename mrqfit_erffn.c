/*************************************************************************** 
MRQFIT_ERFFN: This is the fgauss algorithm taken from Numercial
Recipes but instead of just taking the value of a Gaussian at a given
x value, the Gaussian is intergrated between two x values using the
error function, then divided by the difference between the two x
values to get the mean value. This is usually to take into account the
pixelization of a Gaussian signal. Here's what NR has to say about it:

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
   following equation: y(x) = a0*exp[-(x-a1)^2/2*a2^2] + a3*(x-a1) + a4
5) Have added an additional argument na*2 matrix defining the limits of the
   parameters.
6) If the dyda array is passed as NULL then no attempt is made to fill
   it with the derivatives.
7) Spheres of influence of can be specified for each set of parameters
   (i.e. each Gaussian). If the spheres of influence overlap then for
   the Gaussians the sum is taken, but for the continua an average is
   taken. Sphere of influence for each parameter set should be
   specified in the matrix entry corresponding to the first parameter
   of that set. If you wish to always include a particular parameter
   then either set its bounds to -1 or, if a NULL **inf matrix is
   passed, it is assumed that all sets contribute to all fits. All
   pixels are zero-indexed.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "fit.h"
#include "gamm.h"
#include "memory.h"
#include "const.h"
#include "error.h"

int mrqfit_erffn(double *x, double *a, double **lim, int *ia, int **inf,
		 double *y, double *dyda, int nx, int na, int idx) {

  double bma=0.0,bsmas=0.0,aexp=0.0,bexp=0.0,aexps=0.0,bexps=0.0,expa=0.0,expb=0.0;
  double erfa=0.0,erfb=0.0;
  int    ninc=0,nset=0;
  int    i=0,j=0,k=0;
  int    *inc=NULL;

  /* Check to make sure 2 x limits have been entered */
  if (nx!=2) {
    nferrormsg("mrqfit_erffn(): Only %d x-limits provided.\n\
\tMust provide 2.",nx); return 0;
  }

  /* Check to make sure number of parameters is divisible by 5 */
  if ((MOD(na,5))) {
    nferrormsg("mrqfit_erffn(): Number of input parameters (=%d)\n\
\tnot a multiply of 5",na); return 0;
  }

  /* Check to make sure x-limits are not the same */
  bma=x[1]-x[0]; bsmas=x[1]*x[1]-x[0]*x[0];
  if (x[0]==x[1]) {
    nferrormsg("mrqfit_erffn(): x-limits are equal (=%lf)",x[0]); return 0;
  }

  /** Make some checks on the spheres of influence of each parameter set **/
  /* Initialise inc */
  ninc=0; nset=na/5;
  if ((inc=iarray(nset))==NULL)
    errormsg("mrqfit_erffn(): Cannot allocate memory for inc\n\
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

  /* Add up contributions from all Gaussians */
  for (i=0,k=0,*y=0.0; i<na-1; i+=5,k++) {
    /* Check to make sure width of Gaussian isn't zero */
    if (fabs(a[i+2])<FPMIN) {
      nferrormsg("mrqfit_erffn(): Gaussian %d has near zero\n\
\tor negative width (=%lf)",i/5+1,a[i+2]); return 0;
    }
    /* Force parameters to lie within limits */
    for (j=i; j<i+5; j++) {
      if (ia[j]) { a[j]=(MAX(a[j],lim[j][0])); a[j]=(MIN(a[j],lim[j][1])); }
    }
    /* Calculate function and it's derivatives */
    aexp=(x[0]-a[i+1])/(C_SQRT2*a[i+2]); bexp=(x[1]-a[i+1])/(C_SQRT2*a[i+2]);
    aexps=aexp*aexp; bexps=bexp*bexp;
    expa=(aexps<80.0) ? exp(-aexps) : 0.0; 
    expb=(bexps<80.0) ? exp(-bexps) : 0.0;
    erfa=erffn(aexp); erfb=erffn(bexp);
    *y+=(double)inc[k]*a[i]*a[i+2]*C_SQRTPI/C_SQRT2*(erfb-erfa)/bma + 
      (double)inc[k]/(double)ninc*(a[i+3]/2.0*bsmas/bma-a[i+3]*a[i+1]+a[i+4]);
    if (dyda!=NULL) {
      dyda[i]=a[i+2]*C_SQRTPI/C_SQRT2*(erfb-erfa)/bma;
      dyda[i+1]=a[i]*(expa-expb)/bma-a[i+3];
      dyda[i+2]=a[i]/bma*(C_SQRT2*(aexp*expa-bexp*expb)+
			  C_SQRTPI/C_SQRT2*(erfb-erfa));
      dyda[i+3]=bsmas/2.0/bma-a[i+1];
      dyda[i+4]=1.0;
    }
  }

  free(inc);

  return 1;

}
