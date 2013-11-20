/****************************************************************************
EW(): Calculate the equivalent width in absorption via a choice of
methods. The choice of methods is detailed below. The dispersion can
be zero if the pixels of the provided wavelength, flux, error and
continuum arrays have equal size in wavelength space or disp can be
given in km/s if the pixels have a constant log-linear size. If a NULL
continuum is provided, it is assumed that the supplied flux and error
arrays are normalised. The routine checks for bad pixels in the EW
window. If there are more than nbtol bad pixels, the routine exits
without determining EWs. Pass NBTOL from fit.h to the nbtol parameter
for the default setting. If the verbose flag is set, EW() prints
warnings to stderr, otherwise it is silent. If a status array pixel is
passed to the routine then the pixel status for the EW measurement at
the cwl[0] pixel is returned to the pixel passed (i.e. if cwl[0]
corresponds to pixel i in the spectrum then &(st[i]) should be passed
to EW()).

The status array can currently have the following values: 1 if the EW
measurement is successful and has proceeded as normal. 0 if it was not
possible to make an EW measurement (in this case the EW error is
returned as -1.0 as well). GJSNGLR if during options 2-4 the fit to
the data was degenerate, causing the gaussj() routine inside mrqmin()
to form a singular matrix. MRQITMAX if during options 2-4 the number
of iterations maxed out at MRQNMAX. If the status array is passed as
NULL then this information is not recorded.

opt = 0 : The EW is calculated by simply adding the pixel fluxes
within the specified window (with size win in km/s) with the central
wavelength cwl. Fractional pixels are used at the edges of the window.

opt = 1 : The EW is found by adding the pixel fluxes within the
specified window weighted by a profile function which is taken to be a
pixelized Gaussian with a width given by the instrumental resolving
power, res = lambda/dlambda. Therefore, one must set res to be a
non-zero value when using this option.

opt = 2 : A pixelized Gaussian is fitted to the normalized flux
centered at cwl. The user may specify their own initial guesses for
the fit parameters, as well as their own free parameters and the
limits upon them. If a, lim & ia are passed as NULL then defaults are
used. By default, only the height and the width of the Gaussian are
free parameters. The width of the fitted profile is limited to be
larger than that implied by the instrumental resolving power,
res=lambda/dlambda. Therefore res must be non-zero when using this
option. The maximum width is limited to be less than EW_NWIDTH_FRAC
times the dispersion times the number of pixels in the window passed
to EW (specified by the win parameter). The height of the Gaussian is
limited to be greater than -1.0 (i.e. zero flux). The fit is
terminated when the relative difference in chisq between fitting
iterations falls below EW_CHISQ_PREC. If the user specifies their own
fit parameters then the final fit parameters are passed back out via
the input parameter array, a. The parameter array is passed to
mrqfit_erffn() so a[0] is height, a[1] is position, a[2] is width of
the Gaussian, a[3] gives the slope of the continuum and a[4] gives the
offset of the continuum from zero. If a fit array (of size np) is
passed to the routine then the Gaussian fit is returned, otherwise fit
should be passed as NULL. The returned fit is NOT normalised, unless
the continuum is passed as NULL to EW().

opt = 3 : A pixelised Gaussian doublet is fitted to the normalized
flux, with transitions centred at cwl[0] (lower) and cwl[1]
(upper). The user may specify their own initial guesses for the fit
parameters, as well as their own free parameters and the limits and
constraints upon them. If a, lim & ia are passed as NULL then defaults
are used. By default, only the height and width of each Gaussian are
free parameters, the separation between the two transitions is
constrained and the widths of the two transitions are constrained to
be the same. The default limits on the height and width are the same
as option 2 and the fit is terminated in the same way. The parameter
array is passed to mrqfit_multierffn() so a[0]-a[4] are as specified
for option 2, for the primary (lower wavelength) Gaussian and
a[5]-a[9] are same parameters for the secondary Gaussian. Constraints
are always assumed to be ratios. e.g. if you wish to constrain the
separation of the doublet then initially a[1]=cwl[0], but
a[6]=cwl[1]/cwl[0]. Or if you wish to constrain the widths to be the
same, then a[7]=1.0. For option 3 we assume that the doublet is
unresolved, thus single EW and error is returned for the combined
doublet. If a fit array (of size np) is passed to the routine then the
Gaussian fit is returned, otherwise fit should be passed as NULL. The
returned fit is NOT normalised, unless no continuum is passed to EW().

opt = 4 : This is the same as option 3, except that we assume that the
doublet is resolved, thus separate EWs and errors are returned for
each Gaussian. Also the widths of the two Gaussians are unconstrained
with respect to each other by default. Thus ew and eew must be passed
as a two element array. EW information for the lower wavelength line
is passed back in ew[0] and eew[0] and for the upper wavelength line
it is passed back in ew[1] and eew[1]. If a fit array (of size 2*np) is
passed to the routine then the Gaussian fit for each transition is
returned (in the same array, one after the other), otherwise fit
should be passed as NULL. The returned fit is NOT normalised, unless
no continuum is passed to EW().

opt = 5 : This is like option 3, except that the parameters you pass to
EW() are the final profile shape, i.e. no fit is attempted the
parameters are assumed to be correct.  This is useful when calculating
the sensitivity function, g(z), for a doublet search which has been
conducted using option 3.

opt = 6 : This is like option 4, except that the parameters you pass to
EW() are the final profile shape, i.e. no fit is attempted the
parameters are assumed to be correct.  This is useful when calculating
the sensitivity function, g(z), for a doublet search which has been
conducted using option 4.

cwl contains the first guess at the wavelength of the transition being
fitted. For options 2-4 this guess may be updated by the fitting
routine, if the user wishes (by using user specified fit parameters).
When using options 3-4 cwl must be a two element array containing the
guesses for both the lower (cwl[0]) and upper (cwl[1]) transitions
being fitted.

In all cases the routine returns EW=0 and EEW=-1.0 if there are no
valid pixels in the specified wavelength range. The calling program
must therefore test the EEW returned to determine if the EW is valid at
cwl[0] (& cwl[1]).

Both niter and chisq may be passed as NULL to EW(). For options 2-4, 
if they are not null then the number of iterations and the chi-squared
of the fit are returned to the calling routine.

NOTES 1: This subroutine assumes that the flux is given as F_lambda(lambda)
      2: ia, a & lim should be passed as NULL unless using options 2-4.
      3: If the user wishes to specify their own fit parameters for 
         options 2-4 then the user mush specify all initial fit parameters,
         free parameters, limits and constraints required for that option.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "fit.h"
#include "gamm.h"
#include "memory.h"
#include "utils.h"
#include "const.h"
#include "error.h"

#define EW_INFIN       1.e32
#define EW_NWIDTH_FRAC 1.00
#define EW_CHISQ_PREC  1.e-3

#define FREE_P    free(P);
#define FREE_INFI if (infi!=NULL) { free(*infi); free(infi); }
#define FREE_FIT  if (ia==NULL && a==NULL && lim==NULL) \
                  { free(iai); free(ai); free(*limi); free(limi); }
#define FREE_DATA free(y); free(sig); free(*x); free(x); free(*covar); \
                  free(covar); free(*alpha); free(alpha);
#define FREE_ALL  if (P!=NULL) FREE_P; \
                  if (opt>1) { FREE_FIT; FREE_DATA; FREE_INFI; }

int EW(double *wl, double *fl, double *er, double *co, int np, double *cwl,
       double win, double disp, double res, double *a, double **lim, 
       int *ia, double *ew, double *eew, double *fit, double *chisq,
       int *st, int *niter, int nbtol, int opt, int verbose) {

  double erfden[2],hdisp=0.0,mhdisp=0.0,Psum=0.0;
  double swl[2],ewl[2],dwl=0.0,sp_frac=0.0,ep_frac=0.0;
  double lsp=0.0,lep=0.0,lp=0.0,up=0.0;
  double sum=0.0,esum=0.0,esum_temp=0.0;
  double ochisq=EW_INFIN,chisqi=EW_INFIN,alambda=0.0;
  double *P=NULL;
  double *y=NULL,*sig=NULL,*ai=NULL;
  double **x=NULL,**limi=NULL,**covar=NULL,**alpha=NULL;  
  int    swlidx[2],ewlidx[2],nfit=0;
  int    ret=0; /* Return value from mrqmin */
  int    i=0,j=0;
  int    *iai=NULL,**infi=NULL;

  /* Ensure that at least two pixels are input */
  if (np<2) {
    nferrormsg("EW(): Only %d pixels entered. Need at least 2 pixels.",np);
    return 0;
  }
  /* Make sure window is set with proper width */
  if (win<=0.0) {
    nferrormsg("EW(): Window width is set to %lf. Must be >=0.0",win); 
    return 0;
  }
  /* Check nbtol is sensible */
  if (nbtol<0 || nbtol>np) {
    nferrormsg("EW(): Bad pixel tolerance, nbtol =%d, is nonsensical.",nbtol); 
    return 0;
  }
  /* Check for sensible passing of fit parameter input */
  if (!(ia==NULL && a==NULL && lim==NULL) && 
      !(ia!=NULL && a !=NULL && lim!=NULL))
    errormsg("EW(): Must either pass all or none of\n\
\tia, a and lim as NULL");

  /* Set the +ve and -ve half-dispersions */
  hdisp=(disp) ? 0.5*disp/C_C_K : 0.5*(wl[1]-wl[0]);
  mhdisp=(disp) ? hdisp/(1.0+2.0*hdisp) : hdisp;

  /* Initialise status array if present */
  if (st!=NULL) st[0]=1;

  /* Calculate the denominator for the exponential in the error function */
  if (opt) erfden[0]=C_SQRT2*cwl[0]/(res*C_FWHMSIG);
  else erfden[0]=0.0;
  if (opt==3 || opt==4 || opt==5 || opt==6)
    erfden[1]=C_SQRT2*cwl[1]/(res*C_FWHMSIG);
  else erfden[1]=0.0;

  /** Calculate starting and ending wavelengths for all options **/
  swl[0]=cwl[0]*(1.0-0.5*win/C_C_K); ewl[0]=cwl[0]*(1.0+0.5*win/C_C_K);
  if (opt==3 || opt==4 || opt==5 || opt==6) {
    swl[1]=cwl[1]*(1.0-0.5*win/C_C_K); ewl[1]=cwl[1]*(1.0+0.5*win/C_C_K);
  }
  else { swl[1]=0.0; ewl[1]=0.0; }

  /** Find the corresponding pixels **/
  /* Starting pixel, lower wavelength line. Fraction calculated, but
     only used in option 0 */
  if ((i=idxdval(wl,np,swl[0]))==-1) {
    swlidx[0]=i=np-1; up=(disp) ? wl[i]*(1.0+hdisp) : wl[i]+hdisp;
    if (swl[0]>up) {
      if (verbose) nferrormsg("EW(): Starting wavelength %lf is longer than\n\
\tright edge of final pixel in wl array, %lf",swl[0],up);
      if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0;
      if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; }
      return 1;
    }
    if (!opt) lsp=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
  } else {
    if (!i) {
      swlidx[0]=0; lsp=(disp) ? wl[0]*(1.0-mhdisp) : wl[0]-mhdisp;
      if (swl[0]<lsp) {
	if (verbose) nferrormsg("EW(): Starting wavelength %lf is shorter\n\
\tthan left edge of first pixel in wl array, %lf",swl[0],lsp);
	if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0;
	if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; }
	return 1;
      }
    }
    else {
      if (!opt) {
	lp=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
	if (swl[0]<lp) {
	  swlidx[0]=--i; lsp=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
	}
	else { swlidx[0]=i; lsp=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp; }
      }
      else {
	lp=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
	swlidx[0]=(swl[0]<lp) ? --i : i;
      }
    }
    if (!opt) up=(disp) ? wl[swlidx[0]]*(1.0+hdisp) : wl[swlidx[0]]+hdisp;
  }
  if (!opt) sp_frac=(up-swl[0])/(up-lsp);
  /* Ending pixel, lower wavelength line. Fraction calculated, but
     only used in option 0*/
  if ((i=idxdval(&(wl[swlidx[0]]),np-swlidx[0],ewl[0]))==-1) {
    ewlidx[0]=i=np-1; up=(disp) ? wl[i]*(1.0+hdisp) : wl[i]+hdisp;
    if (ewl[0]>up) {
	if (verbose) nferrormsg("EW(): Ending wavelength %lf is longer than\n\
\tright edge of final pixel in wl array, %lf",ewl[0],up); 
	if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0;
	if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; }
	return 1;
    }
    if (!opt) lep=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp; 
  }
  else {
    i+=swlidx[0];
    if (!i) {
      ewlidx[0]=0; lep=(disp) ? wl[0]*(1.0-mhdisp) : wl[0]-mhdisp;
      if (ewl[0]<lep) {
	if (verbose) nferrormsg("EW(): Ending wavelength %lf is smaller than\n\
\tleft edge of first pixel in wl array, %lf",ewl[0],lep); 
	if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0;
	if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; }
	return 1;
      }
    }
    else {
      if (!opt) {
	lp=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
	if (ewl[0]<lp) {
	  ewlidx[0]=--i; lep=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
	}
	else { ewlidx[0]=i; lep=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp; }
      }
      else {
	lp=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
	ewlidx[0]=(ewl[0]<lp) ? --i : i;
      }
    }
    if (!opt) up=(disp) ? wl[ewlidx[0]]*(1.0+hdisp) : wl[ewlidx[0]]+hdisp;
  }
  if (!opt) {
    if (swlidx[0]==ewlidx[0]) sp_frac=ep_frac=(ewl[0]-swl[0])/(up-lep);
    else ep_frac=(ewl[0]-lep)/(up-lep);
  }
  if (opt==3 || opt==4 || opt==5 || opt==6) {
    /* Starting pixel, higher wavelength line. */
    if ((i=idxdval(wl,np,swl[1]))==-1) {
      swlidx[1]=i=np-1;
      up=(disp) ? wl[i]*(1.0+hdisp) : wl[i]+hdisp;
      if (swl[1]<=up) {
	if (verbose) nferrormsg("EW(): Starting wavelength %lf for second\n\
\ttransition is longer than right edge of final pixel\n\
\tin wl array, %lf",swl[1],up);
	if (st!=NULL) st[0]=0; ew[0]=0.0; eew[0]=-1.0;
	if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; } 
	return 1;
      }      
    } else {
      if (!i) {
	swlidx[1]=0; lsp=(disp) ? wl[0]*(1.0-mhdisp) : wl[0]-mhdisp;
	if (swl[1]<lsp) {
	  if (verbose) nferrormsg("EW(): Starting wavelength %lf for second\n\
\ttransition is shorter than left edge of first pixel\n\
\tin wl array, %lf",swl[1],lsp);
	  if (st!=NULL) st[0]=0; ew[0]=0.0; eew[0]=-1.0;
	  if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; }
	  return 1;
	}
      }
      else {
	lp=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
	swlidx[1]=(swl[1]<lp) ? --i : i;
      }
    }
    /* Ending pixel, higher wavelength line. */
    if ((i=idxdval(&(wl[swlidx[1]]),np-swlidx[1],ewl[1]))==-1) {
      ewlidx[1]=i=np-1; up=(disp) ? wl[i]*(1.0+hdisp) : wl[i]+hdisp;
      if (ewl[1]>up) {
	if (verbose) nferrormsg("EW(): Ending wavelength %lf for second\n\
\ttransition is longer than right edge of final\n\
\tpixel in wl array, %lf",ewl[1],up); 
	if (st!=NULL) st[0]=0; ew[0]=0.0; eew[0]=-1.0;
	if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; }
	return 1;
      }
    }
    else {
      i+=swlidx[1];
      if (!i) {
	ewlidx[1]=0; lep=(disp) ? wl[0]*(1.0-mhdisp) : wl[0]-mhdisp;
	if (ewl[1]<lep) {
	  if (verbose) nferrormsg("EW(): Ending wavelength %lf for second\n\
\ttransition is smaller than left edge of first\n\
\tpixel in wl array, %lf",ewl[1],lep);
	  if (st!=NULL) st[0]=0; ew[0]=0.0; eew[0]=-1.0;
	  if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; }
	  return 1;
	}
      }
      else {
	lp=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
	ewlidx[1]=(ewl[1]<lp) ? --i : i;
      }
    }
  }
  else { swlidx[1]=-1; ewlidx[1]=-1; }

  /** Check the number of bad pixels in the EW window. Exit without
      calculating the EW if there are more than nbtol. **/
  for (i=swlidx[0],nfit=0; i<=ewlidx[0]; i++) 
    if (er[i]>0.0 && ((co!=NULL && co[i]!=0.0) || co==NULL)) nfit++;
  if (nfit<ewlidx[0]-swlidx[0]+1-nbtol) {
    if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0;
    if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; } 
    return 1; 
  }
  if (opt==3 || opt==4 || opt==5 || opt==6) {
    /* Note that if swlidx[1] and ewlidx[0] overlap (e.g. opt==3) some
       pixels are counted twice. This is correct for bad pixel
       determination. */
    for (i=swlidx[1],nfit=0; i<=ewlidx[1]; i++) 
      if (er[i]>0.0 && ((co!=NULL && co[i]!=0.0) || co==NULL)) nfit++;
    if (nfit<ewlidx[1]-swlidx[1]+1-nbtol) {      
      if (st!=NULL) st[0]=0; ew[0]=0.0; eew[0]=-1.0;
      if (opt==4 || opt==6) { ew[1]=0.0; eew[1]=-1.0; }
      return 1;
    }
    /* Now determine nfit correctly for options 3 & 4 */
    for (i=swlidx[0],nfit=0; i<=ewlidx[1]; i++) {
      if ((i<=ewlidx[0] || i>=swlidx[1]) && er[i]>0.0 &&
	  ((co!=NULL && co[i]!=0.0) || co==NULL)) nfit++;
    }
  }

  /* Allocate memory for weighting array for options 1, 2, 3 and 4 */
  if (opt)
    if ((P=darray(np))==NULL) errormsg("EW(): Cannot allocate memory for\n\
\tP array of size %d",np);
  
  /** First option is to calculate the EW just based on a simple sum
      of pixel fluxes across the given EW window. Fractional pixels
      are used at the edges of the window.**/
  if (!opt) {
    /* Do sums to determine EW and its error */
    /* First pixel */
    dwl=(disp) ? wl[swlidx[0]]*(1.0+hdisp)-lsp : wl[swlidx[0]]+hdisp-lsp;
    if (co!=NULL) {
      sum=(er[swlidx[0]]<=0.0 || co[swlidx[0]]==0.0) ? 0.0 :
	(1.0-fl[swlidx[0]]/co[swlidx[0]])*dwl*sp_frac;
      esum_temp=(er[swlidx[0]]<=0.0 || co[swlidx[0]]==0.0) ? 0.0 :
        er[swlidx[0]]/co[swlidx[0]]*dwl*sp_frac;
      esum=esum_temp*esum_temp;
    }
    else {
      sum=(er[swlidx[0]]<=0.0) ? 0.0 : (1.0-fl[swlidx[0]])*dwl*sp_frac;
      esum_temp=(er[swlidx[0]]<=0.0) ? 0.0 : er[swlidx[0]]*dwl*sp_frac;
      esum=esum_temp*esum_temp;
    }
    /* Last pixel */
    if (swlidx[0]!=ewlidx[0]) {
      dwl=(disp) ? wl[ewlidx[0]]*(1.0+hdisp)-lep : wl[ewlidx[0]]+hdisp-lep;
      if (co!=NULL) {
	sum+=(er[ewlidx[0]]<=0.0 || co[ewlidx[0]]==0.0) ? 0.0 :
	  (1.0-fl[ewlidx[0]]/co[ewlidx[0]])*dwl*ep_frac;
	esum_temp=(er[ewlidx[0]]<=0.0 || co[ewlidx[0]]==0.0) ? 0.0 :
	  er[ewlidx[0]]/co[ewlidx[0]]*dwl*ep_frac;
	esum+=esum_temp*esum_temp;
      }
      else {
	sum+=(er[ewlidx[0]]<=0.0) ? 0.0 : (1.0-fl[ewlidx[0]])*dwl*ep_frac;
	esum_temp=(er[ewlidx[0]]<=0.0) ? 0.0 : er[ewlidx[0]]*dwl*ep_frac;
	esum+=esum_temp*esum_temp;
      }
    }
    /* Other pixels */
    for (i=swlidx[0]+1; i<ewlidx[0]; i++) {
      dwl=(disp) ? wl[i]*(hdisp+mhdisp) : hdisp+mhdisp;
      if (co!=NULL) {
	sum+=(er[i]<=0.0 || co[i]==0.0) ? 0.0 : (1.0-fl[i]/co[i])*dwl;
	esum_temp=(er[i]<=0.0 || co[i]==0.0) ? 0.0 : er[i]/co[i]*dwl;
	esum+=esum_temp*esum_temp;
      }
      else {
	sum+=(er[i]<=0.0) ? 0.0 : (1.0-fl[i])*dwl;
	esum_temp=(er[i]<=0.0) ? 0.0 : er[i]*dwl;
	esum+=esum_temp*esum_temp;
      }
    }      
    /* Final EW, error and SNR */
    if (esum>0.0) { *ew=sum; *eew=sqrt(esum); }
    else { if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0; }
  }
  /** The second option is to use a Gaussian with width equal to the
      instrumental resolution as a weighting function for the
      equivalent width and error **/
  else if (opt==1) {
    /* Set the weight profile, find the normalization and normalize
       all weights */
    for (i=swlidx[0]; i<=ewlidx[0]; i++) {
      if (er[i]<=0.0 || (co!=NULL && co[i]==0.0)) P[i]=0.0;      
      else {
	lep=(disp) ? wl[i]*(1.0-mhdisp) : wl[i]-mhdisp;
	up=(disp) ? wl[i]*(1.0+hdisp) : wl[i]+hdisp;
	P[i]=erffn((up-cwl[0])/erfden[0])-erffn((lep-cwl[0])/erfden[0]);
      }
    }
    for (i=swlidx[0],Psum=0.0; i<=ewlidx[0]; i++) Psum+=P[i];
    if (Psum==0.0) {
      if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0;
      FREE_P; return 1;
    }
    for (i=swlidx[0]; i<=ewlidx[0]; i++) P[i]/=Psum;    
    /* Do sums to determine EW and its error */
    for (i=swlidx[0]; i<=ewlidx[0]; i++) {
      dwl=(disp) ? wl[i]*(hdisp+mhdisp) : hdisp+mhdisp;
      if (co!=NULL) {
	sum+=(P[i]<=0.0) ? 0.0 : P[i]*(1.0-fl[i]/co[i])*dwl;
	esum_temp=(P[i]<=0.0) ? 0.0 : P[i]*er[i]/co[i]*dwl;
	esum+=esum_temp*esum_temp; Psum+=P[i]*P[i];
      }
      else {
	sum+=(P[i]<=0.0) ? 0.0 : P[i]*(1.0-fl[i])*dwl;
	esum_temp=(P[i]<=0.0) ? 0.0 : P[i]*er[i]*dwl;
	esum+=esum_temp*esum_temp; Psum+=P[i]*P[i];
      }
    }      
    /* Final EW, error and SNR */
    if (esum>0.0) { *ew=sum/Psum; *eew=sqrt(esum)/Psum; }
    else { if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0; }
  }
  /** The third option is to fit a Gaussian (a pixelised Gaussian
      constructed using the error function) to the flux to use as a
      weighting function for the profile **/  
  else if (opt==2) {
    /* Allocate memory to fitting matrices and arrays etc. */
    if ((y=darray(nfit))==NULL)
      errormsg("EW(): Cannot allocate memory to y array of size %d",nfit);
    if ((sig=darray(nfit))==NULL)
      errormsg("EW(): Cannot allocate memory to sig array of size %d",nfit);
    if ((x=dmatrix(np,2))==NULL)
      errormsg("EW(): Cannot allocate memory to x matrix of size %dx%d",np,2);
    if ((covar=dmatrix(5,5))==NULL)
      errormsg("EW(): Cannot allocate memory to covar matrix of\n\
\tsize %dx%d",5,5);
    if ((alpha=dmatrix(5,5))==NULL)
      errormsg("EW(): Cannot allocate memory to alpha matrix of\n\
\tsize %dx%d",5,5);
    if (ia==NULL && a==NULL && lim==NULL) {
      if ((iai=iarray(5))==NULL)
	errormsg("EW() Cannot allocate memory to iai array of size %d",5);
      if ((ai=darray(5))==NULL)
	errormsg("EW(): Cannot allocate memory to ai array of size %d",5);
      if ((limi=dmatrix(5,2))==NULL)
	errormsg("EW(): Cannot allocate memory to limi matrix of\n\
\tsize %dx%d",5,2);
    }
    /* Fill arrays */
    i=j=0; while (i<ewlidx[0]-swlidx[0]+1) {
      if (er[i+swlidx[0]]>0.0 && 
	  ((co!=NULL && co[i+swlidx[0]]!=0.0) || co==NULL)) {
	x[j][0]=(disp) ? wl[i+swlidx[0]]*(1.0-mhdisp) : 
	  wl[i+swlidx[0]]-mhdisp;
	x[j][1]=(disp) ? wl[i+swlidx[0]]*(1.0+hdisp) : wl[i+swlidx[0]]+hdisp;
	y[j]=(co!=NULL) ? fl[i+swlidx[0]]/co[i+swlidx[0]] : fl[i+swlidx[0]]; 
	sig[j++]=(co!=NULL) ? 
	  er[i+swlidx[0]]/co[i+swlidx[0]] : er[i+swlidx[0]];
      }
      i++;
    }
    /* Set initial fit parameters and limits, if not specified externally */
    if (ia==NULL && a==NULL && lim==NULL) {
      /* Only vary the height and width of Gaussian */
      iai[0]=iai[2]=1; iai[1]=iai[3]=iai[4]=0;
      /* Make initial guesses for parameters and place appropriate
	 limits */
      ai[0]=-0.01; ai[1]=cwl[0]; ai[2]=erfden[0]/C_SQRT2; ai[3]=0.0; ai[4]=1.0;
      limi[0][0]=-1.0; limi[0][1]=EW_INFIN; limi[2][0]=ai[2];
      limi[2][1]=(disp) ? 
	EW_NWIDTH_FRAC*(double)(ewlidx[0]-swlidx[0]+1)*hdisp*cwl[0] :	
	EW_NWIDTH_FRAC*(double)(ewlidx[0]-swlidx[0]+1)*hdisp;
    }
    else {
      iai=ia; ai=a; limi=lim;
    }
    /* Do fit */
    alambda=-1.0;
    for (i=0; i<MRQNMAX; i++) {
      if (!(ret=mrqmin(x,y,sig,nfit,2,ai,limi,iai,NULL,5,covar,alpha,&chisqi,
		  &alambda,mrqfit_erffn))) {
	nferrormsg("EW(): Error fitting pixelised Gaussian to\n\
\tnormalized data");
	return 0;
      }
      if (ret==GJSNGLR && st!=NULL) st[0]=GJSNGLR; 
      if (chisqi<ochisq && (ochisq-chisqi)/(chisqi)<EW_CHISQ_PREC) break;
      ochisq=chisqi;
    }
    if (i==MRQNMAX) {
      if (verbose) warnmsg("EW(): More than %d iterations required to fit\n\
\tpixelised Gaussian at cwl=%lf",MRQNMAX,cwl[0]);
      if (st!=NULL) st[0]=MRQITMAX;
  }
    if (niter!=NULL) *niter=i;
    alambda=0.0;
    if (!(ret=mrqmin(x,y,sig,nfit,2,ai,limi,iai,NULL,5,covar,alpha,&chisqi,
		&alambda,mrqfit_erffn))) {
      nferrormsg("EW(): Error fitting pixelised Gaussian to normalized data\n\
\tfor last time to detemine covariance matrix"); return 0;
    }
    if (ret==GJSNGLR && st!=NULL) st[0]=GJSNGLR; 
    /* Fill fit array and fill it with reconstructed pixelised
       Gaussian */
    if (fit!=NULL) for(i=0;i<np;i++) fit[i]=(co!=NULL) ? co[i] : 1.0;
    for (i=0,Psum=0.0; i<ewlidx[0]-swlidx[0]+1; i++) {
      x[i][0]=(disp) ?
	wl[i+swlidx[0]]*(1.0-mhdisp) : wl[i+swlidx[0]]-mhdisp;
      x[i][1]=(disp) ?
	wl[i+swlidx[0]]*(1.0+hdisp) : wl[i+swlidx[0]]+hdisp;
      if (!mrqfit_erffn(x[i],ai,limi,iai,NULL,&(P[i]),NULL,2,5,i)) {
	nferrormsg("EW(): Cannot reconstruct pixelised Gaussian\n\
\tfrom paramters");
	return 0;
      }
      if (fit!=NULL) fit[i+swlidx[0]]=(co!=NULL) ? P[i]*co[i+swlidx[0]] : P[i];
      if (er[i+swlidx[0]]>0.0 && 
	  ((co!=NULL && co[i+swlidx[0]]!=0.0) || co==NULL))
	{ P[i]=1.0-P[i]; Psum+=P[i]; }
      else P[i]=0.0;
    }
    /* Produce normalized profile and calculate square sum */
    if (Psum==0.0) { 
      if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0; 
      FREE_ALL; return 1; 
    }
    else for (i=0; i<ewlidx[0]-swlidx[0]+1; i++) P[i]/=Psum;
 
    /* Do sums to determine EW and its error */
    for (i=0,sum=0.0,Psum=0.0; i<ewlidx[0]-swlidx[0]+1; i++) {
      dwl=(disp) ? wl[i+swlidx[0]]*(hdisp+mhdisp) : hdisp+mhdisp;
      sum+=(co!=NULL) ? P[i]*(1.0-fl[i+swlidx[0]]/co[i+swlidx[0]])*dwl : 
	P[i]*(1.0-fl[i+swlidx[0]])*dwl; 
      esum_temp=(co!=NULL) ? P[i]*er[i+swlidx[0]]/co[i+swlidx[0]]*dwl :
	P[i]*er[i+swlidx[0]]*dwl;
      esum+=esum_temp*esum_temp; Psum+=P[i]*P[i];
    }    
    /* Final EW, error and chisq */
    *ew=sum/Psum; *eew=sqrt(esum)/Psum;
    if (chisq!=NULL) *chisq=chisqi;
  }
  /* Fit pixelised Gaussian doublet profile for the EW weighting
     function option 3 gives a combined EW measurement for each
     doublet option 4 gives independent EW measurements for the
     doublet */
  else if (opt==3 || opt==4) {   
    /* Allocate memory to fitting matrices and arrays etc. */
    if ((y=darray(nfit))==NULL)
      errormsg("EW(): Cannot allocate memory to y array of size %d",nfit);
    if ((sig=darray(nfit))==NULL)
      errormsg("EW(): Cannot allocate memory to sig array of size %d",nfit);
    if ((x=dmatrix(np,2))==NULL)
      errormsg("EW(): Cannot allocate memory to x matrix of size %dx%d",np,2);
    if ((covar=dmatrix(10,10))==NULL)
      errormsg("EW(): Cannot allocate memory to covar matrix of\n\
\tsize %dx%d",10,10);
    if ((alpha=dmatrix(10,10))==NULL)
      errormsg("EW(): Cannot allocate memory to alpha matrix of\n\
\tsize %dx%d",10,10);
    if ((infi=imatrix(10,2))==NULL)
      errormsg("EW(): Cannot allocate memory to infi matrix of\n\
\tsize %dx%d",10,2);
    if (ia==NULL && a==NULL && lim==NULL) {
      if ((iai=iarray(10))==NULL)
	errormsg("EW() Cannot allocate memory to iai array of size %d",10);
      if ((ai=darray(10))==NULL)
	errormsg("EW(): Cannot allocate memory to ai array of size %d",10);
      if ((limi=dmatrix(10,3))==NULL)
	errormsg("EW(): Cannot allocate memory to limi matrix of\n\
\tsize %dx%d",10,3);
    }        
    /* Fill arrays */
    i=swlidx[0]; j=0; while (i<=ewlidx[1]) {
      if ((i<=ewlidx[0] || i>=swlidx[1]) && 
	  er[i]>0.0 && 
	  ((co!=NULL && co[i]!=0.0) || co==NULL)) {	    
	x[j][0]=(disp) ? wl[i]*(1.0-mhdisp) : 
	  wl[i]-mhdisp;
	x[j][1]=(disp) ? wl[i]*(1.0+hdisp) : wl[i]+hdisp;
	y[j]=(co!=NULL) ? fl[i]/co[i] : fl[i]; 
	sig[j++]=(co!=NULL) ?
	  er[i]/co[i] : er[i];
      }
      i++;
    } 

    /* Set initial fit parameters and limits, if not specified externally */
    if (ia==NULL && a==NULL && lim==NULL) {
      /* Vary only the height and width of the Gaussians */
      iai[0]=iai[2]=1; iai[1]=iai[3]=iai[4]=0;
      iai[5]=iai[7]=1; iai[6]=iai[8]=iai[9]=0;
      for (i=0;i<10;i++) limi[i][2]=0;
      limi[6][2]=1; /* Constrain separation */
      /* Set initial guesses for parameters and limits */
      ai[0]=-0.01; ai[1]=cwl[0]; ai[2]=erfden[0]/C_SQRT2; ai[3]=0.0; ai[4]=1.0;
      ai[5]=-0.01; ai[6]=cwl[1]/cwl[0]; ai[7]=erfden[1]/C_SQRT2; ai[8]=0.0; 
      ai[9]=1.0;
      /* If unres. (opt 3) constrain width of 2ndary transition */
      if (opt==3) {
	limi[7][2]=1;
	iai[7]=0; 
	ai[7]=1.0;
	limi[9][2]=1;	
      }
      limi[0][0]=-1.0; limi[0][1]=EW_INFIN; limi[2][0]=erfden[0]/C_SQRT2;
      limi[2][1]=(disp) ? 
	EW_NWIDTH_FRAC*(double)(ewlidx[0]-swlidx[0]+1)*hdisp*cwl[0] :	
	EW_NWIDTH_FRAC*(double)(ewlidx[0]-swlidx[0]+1)*hdisp;
      limi[5][0]=-1.0; limi[5][1]=EW_INFIN; 
      if (opt==4 ) {
	limi[7][0]=erfden[1]/C_SQRT2;
	limi[7][1]=(disp) ? 
	  EW_NWIDTH_FRAC*(double)(ewlidx[1]-swlidx[1]+1)*hdisp*cwl[1] :	
	  EW_NWIDTH_FRAC*(double)(ewlidx[1]-swlidx[1]+1)*hdisp;
      }
    } else { 
      iai=ia; ai=a; limi=lim; 
    }
    /* Set spheres of influence */
    if (opt==3) {
      infi[0][0]=-1; infi[0][1]=-1; infi[5][0]=-1; infi[5][1]=-1;
    }
    else {
      infi[0][0]=0; infi[0][1]=ewlidx[0]-swlidx[0];
      if (swlidx[1]>ewlidx[0]) {
	infi[5][0]=infi[0][1]+1; infi[5][1]=ewlidx[1]-infi[0][1];
      }
      else {
	infi[5][0]=swlidx[1]-swlidx[0]; infi[5][1]=ewlidx[1]-swlidx[0];
      }
    }     
    /* Do fit */
    alambda=-1.0;
    for (i=0; i<MRQNMAX; i++) {
      if (!(ret=mrqmin(x,y,sig,nfit,2,ai,limi,iai,infi,10,covar,alpha,&chisqi,
		  &alambda,mrqfit_multierffn))) {
	nferrormsg("EW(): Error fitting pixelised Gaussian doublet to\n\
\tnormalized data");
	return 0;
      }
      if (ret==GJSNGLR && st!=NULL) st[0]=GJSNGLR; 
      if (chisqi<ochisq && (ochisq-chisqi)/(chisqi)<EW_CHISQ_PREC) break;
      ochisq=chisqi;
    }
    if (i==MRQNMAX) {
      if (verbose) warnmsg("EW(): More than %d iterations required to fit\n\
\tpixelised Gaussian at cwl=%lf",MRQNMAX,cwl[0]);
      if (st!=NULL) st[0]=MRQITMAX;
    }
    if (niter!=NULL) *niter=i;
    alambda=0.0;
    if (!(ret=mrqmin(x,y,sig,nfit,2,ai,limi,iai,infi,10,covar,alpha,&chisqi,
		&alambda,mrqfit_multierffn))) {
      nferrormsg("EW(): Error fitting pixelised Gaussian doublet to\n\
\tnormalized data for last time to detemine covariance matrix"); return 0;
    }
    if (ret==GJSNGLR && st!=NULL) st[0]=GJSNGLR; 
    /* Fill fit array and fill it with reconstructed pixelised
       Gaussians */
    if (fit!=NULL) for (i=0;i<np;i++) fit[i]=(co!=NULL) ? co[i] : 1.0;
    if (opt==3) {
      for (i=0,Psum=0.0; i<ewlidx[1]-swlidx[0]+1; i++) {
	x[i][0]=(disp) ? wl[i+swlidx[0]]*(1.0-mhdisp) : wl[i+swlidx[0]]-mhdisp;
	x[i][1]=(disp) ? wl[i+swlidx[0]]*(1.0+hdisp) : wl[i+swlidx[0]]+hdisp;
	if (!mrqfit_multierffn(x[i],ai,limi,iai,infi,&(P[i]),NULL,2,10,i)) {
	  nferrormsg("EW(): Cannot reconstruct pixelised Gaussian doublet\n\
\tfrom paramters");
	  return 0;
	}
	if (fit!=NULL) fit[i+swlidx[0]]=(co!=NULL) ?
	  P[i]*co[i+swlidx[0]] : P[i];
	if (er[i+swlidx[0]]>0.0 &&
	    ((co!=NULL && co[i+swlidx[0]]!=0.0) || co==NULL))
	  { P[i]=1.0-P[i]; Psum+=P[i]; }
	else P[i]=0.0;
      }
      /* Produce normalized profile and calculate square sum */
      if (Psum==0.0) {
	if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0;
	FREE_ALL; return 1;
      }
      else for (i=0; i<ewlidx[1]-swlidx[0]+1; i++) P[i]/=Psum;
      
      /* Do sums to determine EW and its error */
      for (i=0,sum=0.0,Psum=0.0; i<ewlidx[1]-swlidx[0]+1; i++) {
	dwl=(disp) ? wl[i+swlidx[0]]*(hdisp+mhdisp) : hdisp+mhdisp;
	sum+=(co!=NULL) ? P[i]*(1.0-fl[i+swlidx[0]]/co[i+swlidx[0]])*dwl :
	  P[i]*(1.0-fl[i+swlidx[0]])*dwl; 
	esum_temp=(co!=NULL) ? P[i]*er[i+swlidx[0]]/co[i+swlidx[0]]*dwl :
	  P[i]*er[i+swlidx[0]]*dwl;
	esum+=esum_temp*esum_temp; Psum+=P[i]*P[i];
      }
      /* Final EW, error and SNR */
      *ew=sum/Psum; *eew=sqrt(esum)/Psum;
      if (chisq!=NULL) *chisq=chisqi;
    }
    else {
      /* Option 4 */
      /* Split profile in two individual Gaussians */     
      /* First Gaussian */
      if (fit!=NULL) for (i=np;i<2*np;i++) fit[i]=(co!=NULL) ? co[i-np] : 1.0;
      for (i=0,Psum=0.0; i<ewlidx[0]-swlidx[0]+1; i++) {
	x[i][0]=(disp) ? wl[i+swlidx[0]]*(1.0-mhdisp) : wl[i+swlidx[0]]-mhdisp;
	x[i][1]=(disp) ? wl[i+swlidx[0]]*(1.0+hdisp) : wl[i+swlidx[0]]+hdisp;
	if (!mrqfit_multierffn(x[i],&(ai[0]),&(limi[0]),&(iai[0]),&(infi[0]),
			       &(P[i]),NULL,2,5,i)) {
	  nferrormsg("EW(): Cannot reconstruct pixelised Gaussian doublet\n\
\tfrom parameters");
	  return 0;
	}
	if (fit!=NULL) fit[i+swlidx[0]]=(co!=NULL) ?
	  P[i]*co[i+swlidx[0]] : P[i];
	if (er[i+swlidx[0]]>0.0 && 
	    ((co!=NULL && co[i+swlidx[0]]!=0.0) || co==NULL)) 
	  { P[i]=1.0-P[i]; Psum+=P[i];}
	else P[i]=0.0;
      }
      
      /* Produce normalized profile and calculate square sum */
      if (Psum==0.0) {
	if (st!=NULL) st[0]=0; ew[0]=ew[1]=0.0; eew[0]=eew[1]=-1.0;
	FREE_ALL; return 1;
      }
      else for (i=0; i<ewlidx[0]-swlidx[0]+1; i++) P[i]/=Psum;
      
      /* Do sums to determine EW and its error */
      esum=esum_temp=dwl=0.0;
      for (i=0,sum=0.0,Psum=0.0; i<ewlidx[0]-swlidx[0]+1; i++) {
	dwl=(disp) ? wl[i]*(hdisp+mhdisp) : hdisp+mhdisp;
	sum+=(co!=NULL) ? P[i]*(1.0-fl[i+swlidx[0]]/co[i+swlidx[0]])*dwl :
	  P[i]*(1.0-fl[i+swlidx[0]])*dwl;
	esum_temp=(co!=NULL) ? P[i]*er[i+swlidx[0]]/co[i+swlidx[0]]*dwl :
	  P[i]*er[i+swlidx[0]]*dwl;
	esum+=esum_temp*esum_temp; Psum+=P[i]*P[i];
      }
      /* Final EW, error and SNR */
      ew[0]=sum/Psum; eew[0]=sqrt(esum)/Psum;

      /* Second Gaussian */
      /* Make parameters as a primary for reconstruction */
      if (limi[5][2]) ai[5]=ai[0]*ai[5];
      if (limi[6][2]) ai[6]=ai[1]*ai[6];
      if (limi[7][2]) ai[7]=ai[2]*ai[7];
      if (limi[8][2]) ai[8]=ai[3]*ai[8];
      if (limi[9][2]) ai[9]=ai[4]*ai[9];
      infi[5][0]=0; infi[5][1]=ewlidx[1]-swlidx[1];
      for (i=0,Psum=0.0; i<ewlidx[1]-swlidx[1]+1; i++) {
	x[i][0]=(disp) ? wl[i+swlidx[1]]*(1.0-mhdisp) : wl[i+swlidx[1]]-mhdisp;
	x[i][1]=(disp) ? wl[i+swlidx[1]]*(1.0+hdisp) : wl[i+swlidx[1]]+hdisp;
	if (!mrqfit_multierffn(x[i],&(ai[5]),&(limi[5]),&(iai[5]),&(infi[5]),
			       &(P[i]),NULL,2,5,i)) {
	  nferrormsg("EW(): Cannot reconstruct pixelised Gaussian doublet\n\
\tfrom parameters");
	  return 0;
	}
	if (fit!=NULL) fit[i+np+swlidx[1]]=(co!=NULL) ?
	  P[i]*co[i+swlidx[1]] : P[i];
	if (er[i+swlidx[1]]>0.0 
	    && ((co!=NULL && co[i+swlidx[1]]!=0.0) || co==NULL))
	  { P[i]=1.0-P[i]; Psum+=P[i];}
	else P[i]=0.0;
      }
      /* Restore parameters to original form */
      if (limi[5][2]) ai[5]=ai[5]/ai[0];
      if (limi[7][2]) ai[7]=ai[7]/ai[2];
      if (limi[6][2]) ai[6]=ai[6]/ai[1];
      if (limi[8][2]) ai[8]=ai[8]/ai[3];
      if (limi[9][2]) ai[9]=ai[9]/ai[4];
	
      /* Produce normalized profile and calculate square sum */
      if (Psum==0.0) { 
	if (st!=NULL) st[0]=0;
	ew[0]=ew[1]=0.0; eew[0]=eew[1]=-1.0;
	FREE_ALL; return 1; 
      }
      else for (i=0; i<ewlidx[1]-swlidx[1]+1; i++) P[i]/=Psum;
      
      /* Do sums to determine EW and its error */
      esum=esum_temp=dwl=0.0;
      for (i=0,sum=0.0,Psum=0.0; i<ewlidx[1]-swlidx[1]+1; i++) {
	dwl=(disp) ? wl[i+swlidx[1]]*(hdisp+mhdisp) : hdisp+mhdisp;
	sum+=(co!=NULL) ? P[i]*(1.0-fl[i+swlidx[1]]/co[i+swlidx[1]])*dwl :
	  P[i]*(1.0-fl[i+swlidx[1]])*dwl;
	esum_temp=(co!=NULL) ? P[i]*er[i+swlidx[1]]/co[i+swlidx[1]]*dwl :
	  P[i]*er[i+swlidx[1]]*dwl;
	esum+=esum_temp*esum_temp; Psum+=P[i]*P[i];
      }
      /* Final EW, error and SNR */
      ew[1]=sum/Psum; eew[1]=sqrt(esum)/Psum;
      if (chisq!=NULL) *chisq=chisqi;
    }
  }
  else if (opt==5 || opt==6) {   
    /* Allocate memory to fitting matrices and arrays etc. */
    if ((y=darray(nfit))==NULL)
      errormsg("EW(): Cannot allocate memory to y array of size %d",nfit);
    if ((sig=darray(nfit))==NULL)
      errormsg("EW(): Cannot allocate memory to sig array of size %d",nfit);
    if ((x=dmatrix(np,2))==NULL)
      errormsg("EW(): Cannot allocate memory to x matrix of size %dx%d",np,2);
    if ((covar=dmatrix(10,10))==NULL)
      errormsg("EW(): Cannot allocate memory to covar matrix of\n\
\tsize %dx%d",10,10);
    if ((alpha=dmatrix(10,10))==NULL)
      errormsg("EW(): Cannot allocate memory to alpha matrix of\n\
\tsize %dx%d",10,10);
    if ((infi=imatrix(10,2))==NULL)
      errormsg("EW(): Cannot allocate memory to infi matrix of\n\
\tsize %dx%d",10,2);

    /* Check that sensible options have been input */
    if (ia==NULL && a==NULL && lim==NULL)
      errormsg("EW(): The user must specify ia, a and lim values for \
option %d",opt);
 
    /* Fill arrays */
    i=swlidx[0]; j=0; while (i<=ewlidx[1]) {
      if ((i<=ewlidx[0] || i>=swlidx[1]) && 
	  er[i]>0.0 && 
	  ((co!=NULL && co[i]!=0.0) || co==NULL)) {	    
	x[j][0]=(disp) ? wl[i]*(1.0-mhdisp) : 
	  wl[i]-mhdisp;
	x[j][1]=(disp) ? wl[i]*(1.0+hdisp) : wl[i]+hdisp;
	y[j]=(co!=NULL) ? fl[i]/co[i] : fl[i]; 
	sig[j++]=(co!=NULL) ?
	  er[i]/co[i] : er[i];
      }
      i++;
    } 

    /* Set initial fit parameters and limits */
    iai=ia; ai=a; limi=lim; 
    
    /* Set spheres of influence */
    if (opt==5) {
      infi[0][0]=-1; infi[0][1]=-1; infi[5][0]=-1; infi[5][1]=-1;
    }
    else {
      infi[0][0]=0; infi[0][1]=ewlidx[0]-swlidx[0];
      if (swlidx[1]>ewlidx[0]) {
	infi[5][0]=infi[0][1]+1; infi[5][1]=ewlidx[1]-infi[0][1];
      }
      else {
	infi[5][0]=swlidx[1]-swlidx[0]; infi[5][1]=ewlidx[1]-swlidx[0];
      }
    }     
    /* Fill fit array and fill it with reconstructed pixelised
       Gaussians */
    if (fit!=NULL) for (i=0;i<np;i++) fit[i]=(co!=NULL) ? co[i] : 1.0;
    if (opt==5) {
      for (i=0,Psum=0.0; i<ewlidx[1]-swlidx[0]+1; i++) {
	x[i][0]=(disp) ? wl[i+swlidx[0]]*(1.0-mhdisp) : wl[i+swlidx[0]]-mhdisp;
	x[i][1]=(disp) ? wl[i+swlidx[0]]*(1.0+hdisp) : wl[i+swlidx[0]]+hdisp;
	if (!mrqfit_multierffn(x[i],ai,limi,iai,infi,&(P[i]),NULL,2,10,i)) {
	  nferrormsg("EW(): Cannot reconstruct pixelised Gaussian doublet\n\
\tfrom paramters");
	  return 0;
	}
	if (fit!=NULL) fit[i+swlidx[0]]=(co!=NULL) ?
	  P[i]*co[i+swlidx[0]] : P[i];
	if (er[i+swlidx[0]]>0.0 &&
	    ((co!=NULL && co[i+swlidx[0]]!=0.0) || co==NULL))
	  { P[i]=1.0-P[i]; Psum+=P[i]; }
	else P[i]=0.0;
      }
      /* Produce normalized profile and calculate square sum */
      if (Psum==0.0) {
	if (st!=NULL) st[0]=0; *ew=0.0; *eew=-1.0;
	FREE_ALL; return 1;
      }
      else for (i=0; i<ewlidx[1]-swlidx[0]+1; i++) P[i]/=Psum;
      
      /* Do sums to determine EW and its error */
      for (i=0,sum=0.0,Psum=0.0; i<ewlidx[1]-swlidx[0]+1; i++) {
	dwl=(disp) ? wl[i+swlidx[0]]*(hdisp+mhdisp) : hdisp+mhdisp;
	sum+=(co!=NULL) ? P[i]*(1.0-fl[i+swlidx[0]]/co[i+swlidx[0]])*dwl :
	  P[i]*(1.0-fl[i+swlidx[0]])*dwl; 
	esum_temp=(co!=NULL) ? P[i]*er[i+swlidx[0]]/co[i+swlidx[0]]*dwl :
	  P[i]*er[i+swlidx[0]]*dwl;
	esum+=esum_temp*esum_temp; Psum+=P[i]*P[i];
      }
      /* Final EW, error and SNR */
      *ew=sum/Psum; *eew=sqrt(esum)/Psum;
    }
    else {
      /* Option 6 */
      /* Split profile in two individual Gaussians */     
      /* First Gaussian */
      if (fit!=NULL) for (i=np;i<2*np;i++) fit[i]=(co!=NULL) ? co[i-np] : 1.0;
      for (i=0,Psum=0.0; i<ewlidx[0]-swlidx[0]+1; i++) {
	x[i][0]=(disp) ? wl[i+swlidx[0]]*(1.0-mhdisp) : wl[i+swlidx[0]]-mhdisp;
	x[i][1]=(disp) ? wl[i+swlidx[0]]*(1.0+hdisp) : wl[i+swlidx[0]]+hdisp;
	if (!mrqfit_multierffn(x[i],&(ai[0]),&(limi[0]),&(iai[0]),&(infi[0]),
			       &(P[i]),NULL,2,5,i)) {
	  nferrormsg("EW(): Cannot reconstruct pixelised Gaussian doublet\n\
\tfrom parameters");
	  return 0;
	}
	if (fit!=NULL) fit[i+swlidx[0]]=(co!=NULL) ?
	  P[i]*co[i+swlidx[0]] : P[i];
	if (er[i+swlidx[0]]>0.0 && 
	    ((co!=NULL && co[i+swlidx[0]]!=0.0) || co==NULL)) 
	  { P[i]=1.0-P[i]; Psum+=P[i];}
	else P[i]=0.0;
      }

      /* Produce normalized profile and calculate square sum */
      if (Psum==0.0) {
	if (st!=NULL) st[0]=0; ew[0]=ew[1]=0.0; eew[0]=eew[1]=-1.0;
	FREE_ALL; return 1;
      }
      else for (i=0; i<ewlidx[0]-swlidx[0]+1; i++) P[i]/=Psum;

      /* Do sums to determine EW and its error */
      esum=esum_temp=dwl=0.0;
      for (i=0,sum=0.0,Psum=0.0; i<ewlidx[0]-swlidx[0]+1; i++) {
	dwl=(disp) ? wl[i]*(hdisp+mhdisp) : hdisp+mhdisp;
	sum+=(co!=NULL) ? P[i]*(1.0-fl[i+swlidx[0]]/co[i+swlidx[0]])*dwl :
	  P[i]*(1.0-fl[i+swlidx[0]])*dwl;
	esum_temp=(co!=NULL) ? P[i]*er[i+swlidx[0]]/co[i+swlidx[0]]*dwl :
	  P[i]*er[i+swlidx[0]]*dwl;
	esum+=esum_temp*esum_temp; Psum+=P[i]*P[i];
      }
      /* Final EW, error and SNR */
      ew[0]=sum/Psum; eew[0]=sqrt(esum)/Psum;

      /* Second Gaussian */
      /* Make parameters as a primary for reconstruction */
      if (limi[5][2]) ai[5]=ai[0]*ai[5];
      if (limi[6][2]) ai[6]=ai[1]*ai[6];
      if (limi[7][2]) ai[7]=ai[2]*ai[7];
      if (limi[8][2]) ai[8]=ai[3]*ai[8];
      if (limi[9][2]) ai[9]=ai[4]*ai[9];
      infi[5][0]=0; infi[5][1]=ewlidx[1]-swlidx[1];
      for (i=0,Psum=0.0; i<ewlidx[1]-swlidx[1]+1; i++) {
	x[i][0]=(disp) ? wl[i+swlidx[1]]*(1.0-mhdisp) : wl[i+swlidx[1]]-mhdisp;
	x[i][1]=(disp) ? wl[i+swlidx[1]]*(1.0+hdisp) : wl[i+swlidx[1]]+hdisp;
	if (!mrqfit_multierffn(x[i],&(ai[5]),&(limi[5]),&(iai[5]),&(infi[5]),
			       &(P[i]),NULL,2,5,i)) {
	  nferrormsg("EW(): Cannot reconstruct pixelised Gaussian doublet\n\
\tfrom parameters");
	  return 0;
	}
	if (fit!=NULL) fit[i+np+swlidx[1]]=(co!=NULL) ?
	  P[i]*co[i+swlidx[1]] : P[i];
	if (er[i+swlidx[1]]>0.0 
	    && ((co!=NULL && co[i+swlidx[1]]!=0.0) || co==NULL))
	  { P[i]=1.0-P[i]; Psum+=P[i]; }
	else P[i]=0.0;
      }
      /* Restore parameters to original form */
      if (limi[5][2]) ai[5]=ai[5]/ai[0];
      if (limi[7][2]) ai[7]=ai[7]/ai[2];
      if (limi[6][2]) ai[6]=ai[6]/ai[1];
      if (limi[8][2]) ai[8]=ai[8]/ai[3];
      if (limi[9][2]) ai[9]=ai[9]/ai[4];
	
      /* Produce normalized profile and calculate square sum */
      if (Psum==0.0) { 
	if (st!=NULL) st[0]=0;
	ew[0]=0.0; eew[0]=-1.0; ew[1]=0.0; eew[1]=-1.0;
	FREE_ALL; return 1; 
      }
      else for (i=0; i<ewlidx[1]-swlidx[1]+1; i++) P[i]/=Psum;

      /* Do sums to determine EW and its error */
      esum=esum_temp=dwl=0.0;
      for (i=0,sum=0.0,Psum=0.0; i<ewlidx[1]-swlidx[1]+1; i++) {
	dwl=(disp) ? wl[i+swlidx[1]]*(hdisp+mhdisp) : hdisp+mhdisp;
	sum+=(co!=NULL) ? P[i]*(1.0-fl[i+swlidx[1]]/co[i+swlidx[1]])*dwl :
	  P[i]*(1.0-fl[i+swlidx[1]])*dwl;
	esum_temp=(co!=NULL) ? P[i]*er[i+swlidx[1]]/co[i+swlidx[1]]*dwl :
	  P[i]*er[i+swlidx[1]]*dwl;
	esum+=esum_temp*esum_temp; Psum+=P[i]*P[i];
      }
      /* Final EW, error and SNR */
      ew[1]=sum/Psum; eew[1]=sqrt(esum)/Psum;
    }
  }
  else errormsg("EW(): Specified option no. %d does not exist",opt);

  /* Clean up */
  FREE_ALL;

  return 1;
}
