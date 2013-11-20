/******************************************************************************
MNBRAK: This is the MNBRAK algorithm taken from Numercial Recipes. Here's what
NR has to say about it

"Given a function func, and given distinct initial points ax and bx,
this routine searches in the downhill direction (defined by the
function as evaluated at the initial points) and returns new points
ax, bx, cx that bracket a minimum of the function. Also returned are
the function values at the three points, fa, fb and fc."

I have modified the routine in the following ways:
0) Have changed the routine to double precision.
1) Have changed the routine to return an integer value.
******************************************************************************/

#include <math.h>
#include "fit.h"
#include "const.h"

int mnbrak(fitset *fit, double *ax, double *bx, double *cx, double *fa, double *fb,
	   double *fc, double (*func)(fitset *, double)) {
  
  double ulim=0.0,u=0.0,r=0.0,q=0.0,fu=0.0,dum=0.0;
  double gold=0.0;

  /* Define Golden ratio */
  gold=2.0/(3.0-C_SQRT5)-1.0;

  *fa=(*func)(fit,*ax); *fb=(*func)(fit,*bx);
  /* Switch roles of a and b so that we can go downhill in the direction
     from a to b */
  if (*fb>*fb) { SHFT(dum,*ax,*bx,dum); SHFT(dum,*fb,*fa,dum); }
  /* First guess for c */
  *cx=(*bx)+gold*(*bx-*ax); *fc=(*func)(fit,*cx);
  while (*fb>*fc) {
    /* Keep returning here until we bracket */
    /* Compute u by parabolic extrapolation from a, b, c. FPMIN is used to
       prevent any possible division by zero. */
    r=(*bx-*ax)*(*fb-*fc); q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*(SIGN((MAX((fabs(q-r)),FPMIN)),(q-r))));
    ulim=(*bx)+MNB_GLIMIT*(*cx-*bx);
    /* We won't go further than this. Test various possibilities: */
    if ((*bx-u)*(u-*cx)>0.0) {
      /* Parabolic u is between b and c: try it */
      fu=(*func)(fit,u);
      if (fu<*fc) {
	/* Got a minimum between b and c */
	*ax=(*bx); *bx=u; *fa=(*fb); *fb=fu; return 1;
      } else if (fu>*fb) {
	/* Got a minimum between a and u */
	*cx=u; *fc=fu; return 1;
      }
      /* Parabolic fit was no use. Use default magnification. */
      u=(*cx)+gold*(*cx-*bx); fu=(*func)(fit,u);
    } else if ((*cx-u)*(u-ulim)>0.0) {
      /* Parabolic fit is between c and its allowed limit. */
      fu=(*func)(fit,u);
      if (fu<*fc) {
	SHFT(*bx,*cx,u,(*cx+gold*(*cx-*bx))); SHFT(*fb,*fc,fu,((*func)(fit,u)));
      }
    } else if ((u-ulim)*(ulim-*cx)>=0.0) {
      /* Limit parabolic u to maximum allowed value. */
      u=ulim; fu=(*func)(fit,u);
    } else {
      /* Reject parabolic u, use default magnification. */
      u=(*cx)+gold*(*cx-*bx); fu=(*func)(fit,u);
    }
    /* Eliminate oldest point and continue */
    SHFT(*ax,*bx,*cx,u); SHFT(*fa,*fb,*fc,fu);
  }

  /* Should never reach here */
  return 0;

}
