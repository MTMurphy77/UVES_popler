/******************************************************************************
 * brent(): Minimization/Maximization of functions. Taken from NR (Ch10.2 2nd *
 *          Ed): Given a function f, and given a bracketing triplet of        *
 *          abscissas ax, bx, cx (such that bx is between ax and cx, and f(bx)*
 *          is less than both f(ax) and f(cx)), this routine isolates the     *
 *          minimum to a fractional precision of about tol using Brent's      *
 *          method. The abscissa of the minimum is returned as xmin, and the  *
 *          minimum function value is returned as brent, the reutrned function*
 *          value.                                                            *
 *                                                                            *
 * NOTES: 1) Have changed the routine to double precision.                    *
 *        2) Have changed brent() to return an integer and pass result out as *
 *           soln.                                                            *
 *                                                                            *
 * Last updated 8th February 2008 by Berkeley Zych                            *
 *****************************************************************************/

#include <math.h>
#include "fit.h"
#include "const.h"
#include "error.h"

int brent(fitset *fit, double ax, double bx, double cx, double (*f)(fitset *,double),
	  double tol, double *xmin, double *soln) {

  double a=0.0,b=0.0,d=0.0,etemp=0.0,fu=0.0,fv=0.0,fw=0.0,fx=0.0,p=0.0,q=0.0,r=0.0;
  double tol1=0.0,tol2=0.0,u=0.0,v=0.0,w=0.0,x=0.0,xm=0.0;
  double e=0.0; /* Distance moved on the step before last */
  double gold=0.0;
  int    i=0;

  /* Define Golden ratio */
  gold=1.0-(3.0-C_SQRT5)/2.0;

  /* a and b must be in ascending order, but input abscissas need not be. */
  a=(ax<cx) ? ax : cx; b=(ax>cx) ? ax : cx;
  /* Initialization */
  x=w=v=bx; fw=fv=fx=(*f)(fit,x);
  /* Main program loop */
  for (i=0; i<BRNMAX; i++) {
    xm=0.5*(a+b); tol2=2.0*(tol1=tol*fabs(x)+BR_ACC);
    /* Test for completion */
    if (fabs(x-xm)<=tol2-0.5*(b-a)) {
      *xmin=x; *soln=fx; return 1;
    }
    if (fabs(e)>tol1) {
      /* Construct a trial parabolic fit */
      r=(x-w)*(fx-fv); q=(x-v)*(fx-fw); p=(x-v)*q-(x-w)*r; q=2.0*(q-r);
      if (q>0.0) p=-p;
      q=fabs(p); etemp=e; e=d;
      if (fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
	d=gold*(e=(x>=xm) ? a-x : b-x);
      /* The above conditions determine the acceptability of the
	 parabolic fit. Here we take the golden section step into the
	 larger of the two segments */
      else {
	/* Take the parabolic step */
	d=p/q; u=x+d; if (u-a<tol2 || b-u<tol2) d=(SIGN(tol1,(xm-x)));
      }
    } else d=gold*(e=(x>=xm) ? a-x : b-x);
    u=(fabs(d)>=tol1) ? x+d : x+(SIGN(tol1,d));
    fu=(*f)(fit,u);
    /* This is the one function evaluation per iteration */
    /* Now decide what to do with our function evaluation */
    /* Housekeeping follows */
    if (fu<=fx) {
      if (u>=x) a=x; else b=x; SHFT(v,w,x,u); SHFT(fv,fw,fx,fu);
    } else {
      if (u<x) a=u; else b=u;
      if (fu<=fw || w==x) { v=w; w=u; fv=fw; fw=fu; }
      else if (fu<=fv || v==x || v==w) { v=u; fv=fu; }
    } /* Done with housekeeping, back for another iteration */
  }
  warnmsg("brent(): Maximum number of iterations (=%d) exceeded",BRNMAX);
  *xmin=x;

  /* Never get here */
  return 0;

}
