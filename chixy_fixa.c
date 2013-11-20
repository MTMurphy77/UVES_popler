/******************************************************************************
 * chixy_fixa() : Taken from NR (Ch15.3 2nd Ed) : Captive function of fitexy, *
 *           returns the value of (chisq - offs) for the slope b=tan(bang).   *
 *           Scaled data and offs are communicated via the global variables.  *
 *                                                                            *
 * Notes : 1) Have changed the routine to double precision.                   *
 *         2) This routine assumes a=const. So does not need changing. There  *
 *            aa must be defined beforehand, externally.                      *
 *                                                                            *
 * Last updated 8th February 2008 by Berkeley Zych                            *
 *****************************************************************************/

#include <math.h>
#include "fit.h"

double chixy_fixa (fitset *fit, double bang) {

  double ans=0.0,b=0.0;
  int    i=0.0;

  b=tan(bang);
  for (i=0; i<fit->n; i++) {
    fit->wx[i]=(SQR((b*fit->sx[i])))+(SQR((fit->sy[i])));
    if (fit->s[i]==1 && fit->c[i]==1)
      fit->wx[i]=(fit->wx[i]<1.0/FIT_INFIN) ? FIT_INFIN : 1.0/fit->wx[i];
  }
  for (ans=-fit->p[1],i=0; i<fit->n; i++) {
    if (fit->s[i]==1 && fit->c[i]==1)
      ans+=fit->wx[i]*(SQR((fit->y[i]-fit->p[0]-b*fit->x[i])));
  }

  return ans;

}
