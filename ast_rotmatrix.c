/******************************************************************************
AST_ROTMATRIX: Compute the precession rotation matrix from the standard
epoch J2000.0 to the specified epoch.
******************************************************************************/

#include "astron.h"

int ast_rotmatrix(double epoch, double **p) {

  double t=0.0,a=0.0,b=0.0,c=0.0,ca=0.0,cb=0.0,cc=0.0,sa=0.0,sb=0.0,sc=0.0;

  /* The rotation matrix coefficients are polynomials in time measured in
     Julian centuries from the standard epoch. The coefficients are in
     radians */
  
  t=(ast_epoch2jd(epoch)-C_JD2000)/(100.0*C_JYEAR);
  a = t*(0.6406161+t*(0.0000839+t*0.0000050))*C_RPDEG;
  b = t*(0.6406161+t*(0.0003041+t*0.0000051))*C_RPDEG;
  c = t*(0.5567530-t*(0.0001185+t*0.0000116))*C_RPDEG;

  /* Compute the cosines and sines once for efficiency */
  ca=cos(a); sa=sin(a); cb=cos(b); sb=sin(b); cc=cos(c); sc=sin(c);

  /* Compute the rotation matrix from the sines and cosines */
  p[0][0]=ca*cb*cc-sa*sb; p[1][0]=-sa*cb*cc-ca*sb; p[2][0]=-cb*sc;
  p[0][1]=ca*sb*cc+sa*cb; p[1][1]=-sa*sb*cc+ca*cb; p[2][1]=-sb*sc;
  p[0][2]=ca*sc; p[1][2]=-sa*sc; p[2][2]=cc;

  return 1;

}
