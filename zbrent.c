/******************************************************************************
 * zbrent(): Taken from NR (Ch9.3 2nd Ed): Using Brent's method, find the     *
 *           root of a function func known to lie between x1 and x2. The root *
 *           returned as zbrent will be refined until its accuracy is tol.    *
 *                                                                            *
 * Notes: 1) Have adapted the routine to return an integer value and pass     *
 *           root back via a pointer instead.                                 *
 *        2) Have changed the routine to double precision.                    *
 *                                                                            *
 * Last updated 7th February 2008 by Berkeley Zych                            *
 *****************************************************************************/

#include <math.h>
#include "fit.h"
#include "error.h"

int zbrent(fitset *fit, double (*func)(fitset *, double), double x1, double x2,
	   double tol, double *root) {
  double a=0.0,b=0.0,c=0.0,d=0.0,e=0.0,min1=0.0,min2=0.0;
  double fa=0.0,fb=0.0,fc=0.0,p=0.0,q=0.0,r=0.0,s=0.0,tol1=0.0,xm=0.0;
  int    i=0;
  
  a=x1; b=c=x2;
  fa=(*func)(fit,a); fb=fc=(*func)(fit,b);
  if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) {
    nferrormsg("zbrent(): Root passed to this routine must be bracketed"); return 0;
  }
  for (i=0; i<ZBNMAX; i++) {
    if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
      /* Rename a, b, c and adjust bounding interval d */
      c=a; fc=fa; e=d=b-a;
    }
    if (fabs(fc)<fabs(fb)) {
      a=b; b=c; c=a; fa=fb; fb=fc; fc=fa;
    }
    /* Convergence check */
    tol1=2.0*ZB_EPS*fabs(b)+0.5*tol; xm=0.5*(c-b);
    if (fabs(xm)<=tol1 || fb==0.0) { *root=b; return 1; }
    if (fabs(e)>=tol1 && fabs(fa)>fabs(fb)) {
      /* Attempt inverse quadratic interpolation */
      s=fb/fa;
      if (a==c) { p=2.0*xm*s; q=1.0-s; }
      else {
	q=fa/fc; r=fb/fc; p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      /* Check whether in bounds */
      if (p>0.0) q=-q; p=fabs(q); min1=3.0*xm*q-fabs(tol1*q); min2=fabs(e*q);
      if (2.0*p<((min1<min2) ? min1 : min2)) {
	/* Accept interpolation */
	e=d; d=p/q;
      } else {
	/* Interpolation failed */
	d=xm; e=d;
      }
    } else {
      /* Bounds decreasing too slowly, use bisection */
      d=xm; e=d;
    }
    /* Move last best guess to a */
    a=b; fa=fb;
    /* Evaluate new trial root */
    if (fabs(d)>tol1) b+=d;
    else b+=(SIGN(tol1,xm));
    fb=(*func)(fit,b);
  }
  warnmsg("zbrent(): Maximum number of iterations (=%d) exceeded",ZBNMAX);

  /* Never reach this point */
  return 1;

}
