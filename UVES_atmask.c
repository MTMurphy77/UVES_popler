/****************************************************************************
* Mask out atmospheric/telluric features in all echelle orders in a
* single 1-D spectrum. The masking is done in the observer frame, so
* the heliocentric correction has to be applied to the atmospheric
* mask wavelength regions if it has been previously applied to the
* spectrum.
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"
#include "const.h"

int UVES_atmask(spectrum *spec, cspectrum *cspec, atmask *amsk, params *par) {

  int      i=0,j=0,k=0;
  int      sidx=0,eidx=0;
  double   lwl=0.0;
  double   *svhwl=NULL,*evhwl=NULL;

  /* Allocate memory for vac.-helio. starting and ending wavelength
     arrays for mask */
  if ((svhwl=darray(amsk->nmask))==NULL)
    errormsg("UVES_atmask(): Cannot allocate memory for svhwl\n\
\tarray of size %d for atmospheric mask from file\n\t\%s",amsk->nmask,amsk->atmaskfile);
  if ((evhwl=darray(amsk->nmask))==NULL)
    errormsg("UVES_atmask(): Cannot allocate memory for evhwl\n\
\tarray of size %d for atmospheric mask from file\n\t\%s",amsk->nmask,amsk->atmaskfile);

  /* Transfer atmospheric mask regions to vac. helio. frame by
     applying same heliocentric correction that was applied to the
     spectrum, if any and if known. */
  if (spec->vhel!=0.0) {
    for (j=0; j<amsk->nmask; j++) {
      svhwl[j]=amsk->swl[j]*(1.0+spec->vhel/C_C_K);
      evhwl[j]=amsk->ewl[j]*(1.0+spec->vhel/C_C_K);
    }
  }

  /** If spectrum is to be masked after being combined with others **/
  if (par->atmask==1) {
    /* Loop over orders in spectrum */
    for (i=0; i<spec->nor; i++) {
      /* Test whether each mask region lies within spectral range of order */
      for (j=0; j<amsk->nmask; j++) {
	lwl=(spec->or[i].csidx>0) ? cspec->rwl[spec->or[i].csidx-1] : cspec->flwl;
	if (svhwl[j]<=cspec->rwl[spec->or[i].ceidx] & evhwl[j]>=lwl) {
	  /* Find the pixels to be masked out */
	  if ((sidx=idxdval(&(cspec->rwl[spec->or[i].csidx]),spec->or[i].nrdp,svhwl[j]))==-1)
	    errormsg("UVES_atmask(): Cannot find pixel in order %d\n\
\t(covering %lf to %lf Ang.) of spectrum %d in file\n\t%s\n\
\twith red edge wavelength >= vac. helio.-corrected starting wavelength,\n\
\t%lf (corrected from %lf by %lf km/s), of atmospheric\n\
\tmask region %d from file\n\t%s",i+1,spec->or[i].fvhlwl,
		     spec->or[i].vhrwl[spec->or[i].np-1],spec->id+1,spec->file,
		     svhwl[j],amsk->swl[j],spec->vhel,j+1,amsk->atmaskfile);
	  if ((eidx=idxdval(&(cspec->rwl[spec->or[i].csidx]),spec->or[i].nrdp,evhwl[j]))==-1)
	    eidx=spec->or[i].nrdp-1;
	  if (sidx>eidx)
	    errormsg("UVES_atmask(): Starting index (%d) greater than\n\
\tending index (%d) for masking atmospheric region %d (from %lf to %lf Ang.,\n\
\theliocentric corrected from %lf to %lf Ang. by %lf km/s) of file\n\t%s\n\
\tin order %d of spectrum %d from file\n\t%s",sidx,eidx,j+1,svhwl[j],evhwl[j],
		     amsk->swl[j],amsk->ewl[j],spec->vhel,amsk->atmaskfile,
		     i+1,spec->id+1,spec->file);
	  /* Now that the pixels to be masked have been found, mask them */
	  for (k=sidx; k<=eidx; k++) spec->or[i].rdst[k]=ACLIP;
	}
      }
    }
  }

  /** If spectrum is to be masked before being combined with others **/
  if (par->atmask==2) {
    /* Loop over orders in spectrum */
    for (i=0; i<spec->nor; i++) {
      /* Test whether each mask region lies within spectral range of order */
      for (j=0; j<amsk->nmask; j++) {
	if (svhwl[j]<=spec->or[i].vhrwl[spec->or[i].np-1] &&
	    evhwl[j]>=spec->or[i].fvhlwl) {
	  /* Find the pixels to be masked out */
	  if ((sidx=idxdval(spec->or[i].vhrwl,spec->or[i].np,svhwl[j]))==-1)
	    errormsg("UVES_atmask(): Cannot find pixel in order %d\n\
\t(covering %lf to %lf Ang.) of spectrum %d in file\n\t%s\n\
\twith red edge wavelength >= vac. helio.-corrected starting wavelength,\n\
\t%lf (corrected from %lf by %lf km/s), of atmospheric\n\
\tmask region %d from file\n\t%s",i+1,spec->or[i].fvhlwl,
		     spec->or[i].vhrwl[spec->or[i].np-1],spec->id+1,spec->file,
		     svhwl[j],amsk->swl[j],spec->vhel,j+1,amsk->atmaskfile);
	  if ((eidx=idxdval(spec->or[i].vhrwl,spec->or[i].np,evhwl[j]))==-1)
	    eidx=spec->or[i].np-1;
	  if (sidx>eidx)
	    errormsg("UVES_atmask(): Starting index (%d) greater than\n\
\tending index (%d) for masking atmospheric region %d (from %lf to %lf Ang.,\n\
\theliocentric corrected from %lf to %lf Ang. by %lf km/s) of file\n\t%s\n\
\tin order %d of spectrum %d from file\n\t%s",sidx,eidx,j+1,svhwl[j],evhwl[j],
		     amsk->swl[j],amsk->ewl[j],spec->vhel,amsk->atmaskfile,
		     i+1,spec->id+1,spec->file);
	  /* Now that the pixels to be masked have been found, mask them */
	  for (k=sidx; k<=eidx; k++) spec->or[i].st[k]=ACLIP;
	}
      }
    }
  }

  /* Clean up */
  free(svhwl); free(evhwl);

  return 1;

}
