/****************************************************************************
* Allocate memory for and initialize wavelength arrays for combined
* spectrum. Requirements are that the number of combined pixels and
* cspec->flwl has been set and par->dv or par->dwl has already been
* set. If opt=0, carry out memory allocation and wavelength array
* initialization. If opt=1 then only carry out memory allocation.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int UVES_init_cspec(cspectrum *cspec, params *par, int opt) {

  double cwls=0.0;
  int    i=0;

  /* Allocate memory for combined spectrum data arrays */
  if ((cspec->wl=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\twavelength array of combined spectrum of size %d",cspec->np);
  if ((cspec->rwl=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tright edge wavelength array of combined spectrum of size %d",cspec->np);
  if ((cspec->fl=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tflux array of combined spectrum of size %d",cspec->np);
  if ((cspec->er=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tflux error array of combined spectrum of size %d",cspec->np);
  if ((cspec->ef=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tflux fluctuation array of combined spectrum of size %d",cspec->np);
  if ((cspec->co=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tcontinuum array of combined spectrum of size %d",cspec->np);
  if ((cspec->res=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tresolution array of combined spectrum of size %d",cspec->np);
  if ((cspec->oco=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\told continuum array of combined spectrum of size %d",cspec->np);
  if ((cspec->no=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tnormalized flux array of combined spectrum of size %d",cspec->np);
  if ((cspec->ne=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tnormalized error array of combined spectrum of size %d",cspec->np);
  if ((cspec->nf=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tnormalized fluctuation array of combined spectrum of size %d",cspec->np);
  if ((cspec->csq=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tRMS array of combined spectrum of size %d",cspec->np);
  if ((cspec->ccsq=darray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tclipped RMS array of combined spectrum of size %d",cspec->np);
  if ((cspec->ncb=iarray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tcomb. number array of combined spectrum of size %d",cspec->np);
  if ((cspec->nccb=iarray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tclipped comb. number array of combined spectrum of size %d",cspec->np);
  if ((cspec->st=iarray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\tstatus array of combined spectrum of size %d",cspec->np);
  if ((cspec->ost=iarray(cspec->np))==NULL)
    errormsg("UVES_init_cspec(): Cannot allocate memory for\n\
\told status array of combined spectrum of size %d",cspec->np);

  /* Initialise continuum and status values */
  for (i=0; i<cspec->np; i++) { cspec->co[i]=-INFIN; cspec->st[i]=NCLIP; }

  if (opt) return 1;

  /* Fill wavelength and right edge wavelength arrays */
  if (par->linear) {
    cwls=cspec->flwl+0.5*cspec->dwl;
    for (i=0; i<cspec->np; i++) {
      cspec->wl[i]=cwls+cspec->dwl*(double)i;
      cspec->rwl[i]=cspec->wl[i]+0.5*cspec->dwl;
    }
  }
  else {
    cwls=cspec->flwl*pow(10.0,0.5*cspec->dv);
    for (i=0; i<cspec->np; i++) {
      cspec->wl[i]=cwls*pow(10.0,cspec->dv*(double)i);
      cspec->rwl[i]=cwls*pow(10.0,cspec->dv*(0.5+(double)i));
    }
  }

  return 1;

}
