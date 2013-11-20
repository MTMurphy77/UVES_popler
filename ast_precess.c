/******************************************************************************
AST_PRECESS: Precess coordinates from epoch1 to epoch2.

The method used here is based on the new IAU system described in the
supplement to the 1984 Astronomical Almanac.  The precession is
done in two steps; precess epoch1 to the standard epoch J2000.0 and then
precess from the standard epoch to epoch2.  The precession between
any two dates is done this way because the rotation matrix coefficients
are given relative to the standard epoch.

Right ascensions are in hours and declinations are in degrees.

******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "astron.h"
#include "memory.h"
#include "error.h"

int ast_precess(double ra1, double dec1, double epoch1, double epoch2,
		double *ra2, double *dec2) {

  double r0[3],r1[3],**p=NULL;

  /* Allocate memory to p matrix */
  if ((p=dmatrix(3,3))==NULL)
    errormsg("ast_precess(): Cannot allocate memory to p 3x3 matrix");

  /* If the input epoch is 0 then assume the input epoch is the same as the
     output epoch. If the two epochs are the same then return the
     coordinates from epoch1. */
  if (epoch1==0.0 || epoch1==epoch2) {
    *ra2=ra1; *dec2=dec1; return 1;
  }

  /* Rectangular equitorial coordinates (direction cosines) */
  *ra2=ra1*15.0*C_RPDEG; *dec2=dec1*C_RPDEG;
  r0[0]=cos(*ra2)*cos(*dec2); r0[1]=sin(*ra2)*cos(*dec2); r0[2]=sin(*dec2);

  /* If epoch1 is not the standard epoch then precess to the standard epoch */
  if (epoch1!=C_J2000) {
    if (!ast_rotmatrix(epoch1,p))
      errormsg("Unknown error returned from ast_rotmatrix()");
    /* Multiply by the inverse of p which is also the transpose of p */
    r1[0]=p[0][0]*r0[0]+p[0][1]*r0[1]+p[0][2]*r0[2];
    r1[1]=p[1][0]*r0[0]+p[1][1]*r0[1]+p[1][2]*r0[2];
    r1[2]=p[2][0]*r0[0]+p[2][1]*r0[1]+p[2][2]*r0[2];
    r0[0]=r1[0]; r0[1]=r1[1]; r0[2]=r1[2];
  }

  /* If epoch2 is not the standard epoch then precess from the standard
     epoch to the desired epoch */
  if (epoch2!=C_J2000) {
    if (!ast_rotmatrix(epoch2,p))
      errormsg("Unknown error returned from ast_rotmatrix()");
    /* Multiply by the inverse of p which is also the transpose of p */
    r1[0]=p[0][0]*r0[0]+p[1][0]*r0[1]+p[2][0]*r0[2];
    r1[1]=p[0][1]*r0[0]+p[1][1]*r0[1]+p[2][1]*r0[2];
    r1[2]=p[0][2]*r0[0]+p[1][2]*r0[1]+p[2][2]*r0[2];
    r0[0]=r1[0]; r0[1]=r1[1]; r0[2]=r1[2];
  }

  /* Convert from radians to hours and degrees */
  *ra2=atan2(r0[1],r0[0])/15.0/C_RPDEG; *dec2=asin(r0[2])/C_RPDEG;
  if (*ra2<0.0) *ra2+=24.0;

  /* Clean up */
  free(*p); free(p);

  return 1;

}
