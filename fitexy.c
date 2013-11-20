/******************************************************************************
FITEXY: This is the FITEXY algorithm taken from Numercial Recipes. Here's what
NR has to say about it

"Straight-line fit to input data * x[1..ndat] and y[1..ndat] with
errors in both x and y, the respective standard deviations being the
input quantities sigx[1..ndat] and sigy[1..ndat]. Output quantities
are a and b such that y = a + bx minimizes chisq, whose value is
returned as chi2. Standard errors on a and b are returned as siga and
sigb. These are not meaningful if either (i) the fit is poor, or (ii)
b is so large that the data are consistent with a vertical (infinite
b) line. If siga and sigb are returned as BIG, then the data are
consistent with all values of b."

I have modified the routine in the following ways:
0) Have changed the routine to double precision.
1) Have adapted the routine to return an integer value.
2) Have added options flag to hold some line parameters constant.
3) Added status arrays for the x and y arrays (sts=1 for valid array elements).
4) Have invoked an iterative rejection algorithm whereby the user
   specifies the number of iterations and sigma rejection
   threshold. If no rejection is to be performed then the user should
   input 0 or 1 for nclip.
5) Options flag opt1: (0) Fits just b, keeping a const at input value;
                      (1) Fits both a & b
6) Options flag opt2: (0) Allow routine to determine initial guess for b;
                      (1) Initial guess b must be same sign as input;
                      (2) Input value of b taken as initial guess
7) If at some stage in the iterative rejection process the number of
   unclipped pixels drops below 2, an error message is reported, the
   algorithm stops and the iteration number is returned as a negative
   integer.
8) Added verbosity flag: (0) Quiet; (1) Verbose;
*****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "fit.h"
#include "memory.h"
#include "stats.h"
#include "const.h"
#include "error.h"

#define FREE_FIT free(fit.x); free(fit.y); free(fit.sx); free(fit.sy); free(fit.wx); \
  free(fit.p); free(fit.s); free(fit.c);
#define DEFAULT_RES { if (opt1) *a=0.0; *b=*siga=*sigb=0.0; }
#define TEST_NCLIP if (fit.nc<2) { \
    if (verb) nferrormsg("fitexy(): Less than 2 valid data points in input arrays"); \
    *b/=scale; if (opt1) *a/=scale; FREE_FIT; return -j; \
  }

int fitexy(double *x, double *y, double *sigx, double *sigy, int *stsx, int *stsy,
	   int ndat, double sigclip, int nclip, double *a, double *b, double *siga,
	   double *sigb, double *chi2, int opt1, int opt2, int verb) {

  double   ang[6],ch[6];
  double   amx=0.0,amn=0.0,varx=0.0,vary=0.0,scale=1.0,bin=0.0,bmn=0.0,bmx=0.0;
  double   d1=0.0,d2=0.0,r2=0.0,sc2=0.0,sx2=0.0,sy2=0.0,x1=0.0,y1=0.0,dx=0.0,dy=0.0;
  double   dum1=0.0,dum2=0.0,dum3=0.0,dum4=0.0,tempr=0.0;
  int      ncit=0,ci=0,clip=0;
  int      i=0,j=0;
  statset  stat;
  fitset   fit;

  /* Ensure number of clipping iterations nclip is sensible */
  ncit=(MAX(1,nclip)); sc2=sigclip*sigclip;

  /* Allocate memory for working arrays */
  if ((fit.x=darray(ndat))==NULL)
    errormsg("fitexy(): Cannot allocate memory for fit.x array of size %d",ndat);
  if ((fit.y=darray(ndat))==NULL)
    errormsg("fitexy(): Cannot allocate memory for fit.y array of size %d",ndat);
  if ((fit.sx=darray(ndat))==NULL)
    errormsg("fitexy(): Cannot allocate memory for fit.sx array of size %d",ndat);
  if ((fit.sy=darray(ndat))==NULL)
    errormsg("fitexy(): Cannot allocate memory for fit.sy array of size %d",ndat);
  if ((fit.wx=darray(ndat))==NULL)
    errormsg("fitexy(): Cannot allocate memory for fit.wx array of size %d",ndat);
  if ((fit.p=darray(2))==NULL)
    errormsg("fitexy(): Cannot allocate memory for fit.p array of size %d",2);
  if ((fit.s=iarray(ndat))==NULL)
    errormsg("fitexy(): Cannot allocate memory for fit.s array of size %d",ndat);
  if ((fit.c=iarray(ndat))==NULL)
    errormsg("fitexy(): Cannot allocate memory for fit.c array of size %d",ndat);

  /* Copy valid data to working arrays */
  fit.n=ndat; for (i=0,fit.ns=0,fit.nc=0; i<fit.n; i++) {
    fit.x[i]=x[i]; fit.sx[i]=sigx[i]; fit.y[i]=y[i]; fit.sy[i]=sigy[i]; 
    if ((stsx==NULL || stsx[i]==1) && (stsy==NULL || stsy[i]==1)) {
      fit.s[i]=fit.c[i]=1; fit.ns++; fit.nc++;
    } else fit.s[i]=fit.c[i]=0;
  }
  TEST_NCLIP;

  /* Find the x and y RMSs */
  if (!stats(fit.x,fit.sx,NULL,NULL,fit.s,fit.n,1,&stat)) {
    nferrormsg("fitexy(): Error returned from stats()\n\
\twhen calculating variance of x array"); DEFAULT_RES; FREE_FIT; return 0;
  }
  /* varx=stat.emean; */
  varx=stat.ewmean;
  if (!stats(fit.y,fit.sy,NULL,NULL,fit.s,fit.n,1,&stat)) {
    nferrormsg("fitexy(): Error returned from stats()\n\
\twhen calculating variance of y array"); DEFAULT_RES; FREE_FIT; return 0;
  }
  /* vary=stat.emean; */
  vary=stat.ewmean;
  /* Scale the data */
  scale=varx/vary;
  for (i=0; i<fit.n; i++) {
    if (fit.s[i]==1) {
      fit.y[i]*=scale; fit.sy[i]*=scale;
      /* Use both x and y weights in first trial fit */
      fit.wx[i]=sqrt((SQR(fit.sx[i]))+(SQR(fit.sy[i])));
    }
  }
  /* Trial fit for b */
  bin=*b;
  if (opt2<=1) {
    if (!opt1) dum1=(*a);
    if (!linfit(fit.x,fit.y,fit.wx,fit.s,fit.n,1,opt1,&dum1,b,&dum2,&dum3,&dum4)) {
      nferrormsg("fitexy(): Error returned from linfit()"); DEFAULT_RES; FREE_FIT;
      return 0;
    }
  }
  if (opt2==1 && bin!=0.0 && *b*bin<0.0) *b=bin;
  /* Start sigma clipping iterations with unclipped arrays */
  while (j<ncit) {
    /* Increase iteration counter */
    j++;
    /* Construct several angles for reference points and make b an angle */
    fit.p[1]=ang[0]=0.0; ang[1]=atan(*b); ang[3]=0.0; ang[4]=ang[1]; ang[5]=C_PI_2;
    /* Fully varying a and b */
    if (opt1) {
      for (i=3; i<6; i++) ch[i]=chixy(&fit,ang[i]);
      /* Bracket the chisq minimum and then locate it with brent */
      mnbrak(&fit,&ang[0],&ang[1],&ang[2],&ch[0],&ch[1],&ch[2],chixy);
      if (!brent(&fit,ang[0],ang[1],ang[2],chixy,FXY_ACC,b,chi2)) {
	nferrormsg("fitexy(): Error returned from brent()"); DEFAULT_RES; FREE_FIT;
	return 0;
      }
      *chi2=chixy(&fit,*b); *a=fit.p[0];
    } else {
      /* Fixed a, vary b */
      fit.p[0]=(*a)*scale; for (i=3; i<6; i++) ch[i]=chixy_fixa(&fit,ang[i]);
      /* Bracket the chisq minimum and then locate it with brent */
      mnbrak(&fit,&ang[0],&ang[1],&ang[2],&ch[0],&ch[1],&ch[2],chixy_fixa);
      if (!brent(&fit,ang[0],ang[1],ang[2],chixy_fixa,FXY_ACC,b,chi2)) {
	nferrormsg("fitexy(): Error returned from brent()"); DEFAULT_RES; FREE_FIT;
	return 0;
      }
      *chi2=chixy_fixa(&fit,*b);
    }
    /* Save the inverse sum of weights at the minimum */
    for (i=0,r2=0.0; i<fit.n; i++) if (fit.c[i]==1) r2+=fit.wx[i];
    r2=1.0/r2;
    /* Now find the standard errors for b as points where delta(chisq) = 1 */
    bmx=bmn=FIT_INFIN; fit.p[1]=(*chi2)+1.0;
    for (i=0; i<6; i++) {
      /* Go through saved values to bracket the desired roots.  Note
	 periodicity in slope angles */
      if (ch[i]>fit.p[1]) {
	d1=fabs(ang[i]-(*b)); while (d1>=C_PI) d1-=C_PI; d2=C_PI-d1;
	if (ang[i]<*b) { SWAP(d1,d2); }
	if (d1<bmx) bmx=d1;
	if (d2<bmn) bmn=d2;
      }
    }
    /* Call zbrent to find the roots */
    if (bmx<FIT_INFIN) {
      if (opt1) {
	/* Fully vary a and b */
	if (!zbrent(&fit,chixy,*b,*b+bmx,FXY_ACC,&bmx)) {
	  nferrormsg("fitexy(): Error returned from zbrent()"); DEFAULT_RES; FREE_FIT;
	  return 0;
	}
	bmx-=(*b); amx=fit.p[0]-(*a);
	if (!zbrent(&fit,chixy,*b,*b-bmn,FXY_ACC,&bmn)) {
	  nferrormsg("fitexy(): Error returned from zbrent()"); DEFAULT_RES; FREE_FIT;
	  return 0;
	}
	bmn-=(*b); amn=fit.p[0]-(*a);
	*sigb=sqrt(0.5*((SQR(bmx))+(SQR(bmn))))/(scale*(SQR((cos(*b))))); 
	*siga=sqrt(0.5*((SQR(amx))+(SQR(amn)))+r2)/scale;
      } else {
	/* Fixed a, vary b */
	if (!zbrent(&fit,chixy_fixa,*b,*b+bmx,FXY_ACC,&bmx)) {
	  nferrormsg("fitexy(): Error returned from zbrent()"); DEFAULT_RES; FREE_FIT;
	  return 0;
	}
	bmx-=(*b);
	if (!zbrent(&fit,chixy_fixa,*b,*b-bmn,FXY_ACC,&bmn)) {
	  nferrormsg("fitexy(): Error returned from zbrent()"); DEFAULT_RES; FREE_FIT;
	  return 0;
	}
	bmn-=(*b);
	*sigb=sqrt(0.5*((SQR(bmx))+(SQR(bmn))))/(scale*(SQR((cos(*b))))); 
	*siga=0.0;
      }
    } else (*sigb)=(*siga)=FIT_INFIN;
    /* Convert b from being an angle */
    *b=tan(*b);
    /* Clip array if necessary */
    if (j<ncit) {
      for (i=0,clip=0,fit.nc=0; i<fit.n; i++) {
	if (fit.s[i]==1) {
	  sx2=1.0/fit.sx[i]/fit.sx[i]; sy2=1.0/fit.sy[i]/fit.sy[i];
	  x1=(fit.x[i]*sx2+(*b*(fit.y[i]-(*a)))*sy2)/(sx2+(*b)*(*b)*sy2);
	  y1=*a+(*b)*x1; dx=x1-fit.x[i]; dy=y1-fit.y[i]; ci=fit.c[i];
	  if (dx*dx*sx2+dy*dy*sy2>sc2) fit.c[i]=0;
	  else { fit.c[i]=1; fit.nc++; }
	  if (!clip && fit.c[i]!=ci) clip=1;
	}
      }
      TEST_NCLIP;
      if (!clip) j=ncit;
    }
  }
  /* Unscale the answers */
  *b/=scale; if (opt1) *a/=scale;

  /* Clean up */
  FREE_FIT;

  return 1;

}
