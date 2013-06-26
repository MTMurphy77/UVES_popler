/****************************************************************************
* Calculate, using a simply motivated model from Hans Dekker (designer
* of UVES), the angular CCD pixel scale, as projected on the sky, as a
* function of the input pixel x-coordinate (i.e. CCD dispersion
* coordinate). The input value of cwl (central wavelength of setting)
* is only used to distinguish between the blue and red arms of UVES
* (which have quite different camera focal lengths). The number of
* pixels in the echelle order, n, should be input by the user and it
* is assumed that this number of pixels in centred on the echellogram
* on the CCD. The CCD binning factor in the dispersion direction
* should also input.
*
* The pixel scale is given by
*
* pixscal = 3600.0 * arctan[Spix*Dcoll/Fcam/M/Dtel] where
* Spix = Physical pixel size,
* Dcoll = Collimated beam diameter,
* Fcam = Camera focal length,
* Dtel = Telescope diamter,
* M = Anamorphic magnification factor (depends on position along order).
* M is given by Cos(alpha)/Cos(beta) where alpha is slightly different
* for the two gratings and beta is the angle of incidence (changes
* along order). Values for these different parameters are hardcoded in
* mm and degrees here.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "const.h"

double UVES_pixscal(double cwl, double x, int n, int binx) {

  double Spix=0.015,Dcoll=200.0,Dtel=8000.0;
  double Fcam[2]={360.0,500.0};
  double alpha[2]={76.0,75.04};
  double pixscal=0.0,alp=0.0,bet=0.0,MFcam=0.0;

  if (cwl<UVESBORR) {
    alp=C_RPDEG*alpha[0]; bet=alp+atan(Spix*(double)binx*(x-(double)(n/2))/Fcam[0]);
    MFcam=Fcam[0]*cos(alp)/cos(bet);
  } else {
    alp=C_RPDEG*alpha[1]; bet=alp+atan(Spix*(double)binx*(x-(double)(n/2))/Fcam[1]);
    MFcam=Fcam[1]*cos(alp)/cos(bet);
  }
  pixscal=C_ASPDEG*atan(Spix*Dcoll/MFcam/Dtel)*(double)binx/C_RPDEG;

  return pixscal;

}
