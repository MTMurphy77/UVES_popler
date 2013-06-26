/****************************************************************************
* Combined the continuua from all orders of all spectra into one
* continuous, smooth continuum for the combined spectrum
*
* opt = 0: Scale all orders according to their exposure time
* opt = 1: Do no orders according to their exposure times
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "stats.h"
#include "memory.h"
#include "error.h"

int UVES_combine_cont(spectrum *spec, int nspec, cspectrum *cspec, int opt,
		      params *par) {

  double   scale=0.0,slope=0.0;
  double   *dat=NULL,*wgt=NULL;
  int      ndat=0,maxndat=0;
  int      lidx=0,ridx=0;
  int      i=0,j=0,k=0,l=0;
  statset  stat;

  /* Scale each spectrum+continuum according to its exposure time */
  maxndat=0; for (i=0; i<nspec; i++) {
    scale=3600.0/spec[i].etime;
    for (j=0; j<spec[i].nor; j++) {
      /* Check first to see if order is useful */
      if (spec[i].or[j].nuse>=MINUSE) {
	if (!opt) {
	  for (k=0; k<spec[i].or[j].nrdp; k++) {
	    spec[i].or[j].rdfl[k]*=scale; spec[i].or[j].rdco[k]*=scale;
	    if (spec[i].or[j].rdst[k]==1) {
	      spec[i].or[j].rder[k]*=scale; spec[i].or[j].rdef[k]*=scale;
	      spec[i].or[j].rdme[k]*=scale;
	    }
	  }
	}
	maxndat++;
      }
    }
  }

  /* Find maximum number of possible contributions to continuum at each
     point in combined spectrum and allocate memory to data arrays */
  maxndat*=2;
  if ((dat=darray(maxndat))==NULL)
    errormsg("UVES_combine_cont(): Could not allocate\n\
\tmemory for dat array of size %d",maxndat);
  if ((wgt=darray(maxndat))==NULL)
    errormsg("UVES_combine_cont(): Could not allocate\n\
\tmemory for wgt array of size %d",maxndat);

  /** Run through combined spectrum to produce combined continuum **/
  for (l=0; l<cspec->np; l++) {
    ndat=0;
    /* Find which pixels in each order contribute to continuum at this point */
    for (i=0; i<nspec; i++) {
      for (j=0; j<spec[i].nor; j++) {
	/* Check first to see if order is useful */
	if (spec[i].or[j].nuse>=MINUSE) {
	  if (spec[i].or[j].csidx<=l && spec[i].or[j].ceidx>=l) {
	    dat[ndat]=spec[i].or[j].rdco[l-spec[i].or[j].csidx];
	    lidx=l-spec[i].or[j].csidx; ridx=spec[i].or[j].ceidx-l;
	    if (lidx<par->contwgt && lidx<spec[i].or[j].nrdp/2)
	      wgt[ndat]=sqrt((double)par->contwgt/(double)(lidx+1));
	    else if (ridx<par->contwgt && ridx<=spec[i].or[j].nrdp/2)
	      wgt[ndat]=sqrt((double)par->contwgt/(double)(ridx+1));
	    else wgt[ndat]=1.0;
	    ndat++;
	  }
	}
      }
    }
    /* Assign a mean continuum level for the relevant cspec pixel */
    if (!ndat) cspec->co[l]=-INFIN;
    else if (ndat==1) cspec->co[l]=dat[0]; 
    else {
      /* Compute the weighted mean continuum level */
      if (!stats(dat,wgt,NULL,NULL,NULL,ndat,1,&stat))
	errormsg("UVES_combine_cont(): Unknown error returned\n\
\tfrom stats()");
      cspec->co[l]=stat.wmean;
    }
  }

  /** Run through again, filling in the gaps **/
  l=0; while (l<cspec->np) {
    if (cspec->co[l]==-INFIN) {
      /* Find next valid continuum pixel */
      i=l+1; while (i<cspec->np && cspec->co[i]==-INFIN) i++;
      /* Fill in gap with linear approximation */
      if (i==cspec->np)
	for (j=l; j<i; j++) cspec->co[j]=cspec->co[l-1];
      else {
	slope=(cspec->co[i]-cspec->co[l-1])/(double)(i-l+1);
	for (j=l; j<i; j++) {
	  cspec->co[j]=cspec->co[l-1]+slope*(double)(j-l+1);
	}
      }
      l=i;
    }
    l++;
  }

  /* Scale each spectrum+continuum according to new combined continuum */
  for (i=0; i<nspec; i++) {
    for (j=0; j<spec[i].nor; j++) {
      /* Check first to see if order is useful */
      if (spec[i].or[j].nuse>=MINUSE) {
	for (k=0; k<spec[i].or[j].nrdp; k++) {
	  scale=cspec->co[spec[i].or[j].csidx+k]/spec[i].or[j].rdco[k];
	  spec[i].or[j].rdco[k]=cspec->co[spec[i].or[j].csidx+k];
	  spec[i].or[j].rdfl[k]*=scale; spec[i].or[j].rder[k]*=scale;
	  spec[i].or[j].rdef[k]*=scale;
	}
      }
    }
  }

  /* Clean up */
  free(dat); free(wgt);
  
  return 1;

}
