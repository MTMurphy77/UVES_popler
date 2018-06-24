/****************************************************************************
* Combine specific orders of an array of spectra which have not been
* continuum fitted (but which should have been rescaled to closely match
* each other) into a combined spectrum.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "stats.h"
#include "memory.h"
#include "error.h"

int UVES_combine_nocont(spectrum *spec, int nspec, cspectrum *cspec, int **comb,
			int ncomb, int csidx, int ceidx, params *par) {
  
  double   bigerr=0.0,wgtsq=0.0,sumresonwgtsq=0.0,sum1onwgtsq=0.0;
  double   *dat=NULL,*err=NULL,*efl=NULL,*med=NULL,*res=NULL;
  int      ndat=0,cst=0;
  int      lval=0,rval=0;
  int      post074clip=1;
  int      i=0,j=0,l=0;
  int      *clip=NULL,**con=NULL;
  statset  stat;

  /* Allocate memory to data arrays */
  if ((dat=darray(ncomb))==NULL)
    errormsg("UVES_combine_nocont(): Could not allocate\n\
\tmemory for dat array of size %d",ncomb);
  if ((err=darray(ncomb))==NULL)
    errormsg("UVES_combine_nocont(): Could not allocate\n\
\tmemory for err array of size %d",ncomb);
  if ((efl=darray(ncomb))==NULL)
    errormsg("UVES_combine_nocont(): Could not allocate\n\
\tmemory for efl array of size %d",ncomb);
  if ((med=darray(ncomb))==NULL)
    errormsg("UVES_combine_nocont(): Could not allocate\n\
\tmemory for med array of size %d",ncomb);
  if ((res=darray(ncomb))==NULL)
    errormsg("UVES_combine_nocont(): Could not allocate\n\
\tmemory for res array of size %d",ncomb);
  if ((clip=iarray(ncomb))==NULL)
    errormsg("UVES_combine_nocont(): Could not allocate\n\
\tmemory for clip array of size %d",ncomb);
  if ((con=imatrix(ncomb,3))==NULL)
    errormsg("UVES_combine_nocont(): Could not allocate\n\
\tmemory for con matrix of size %dx%d",ncomb,3);

  /* Loop over all relevant combined pixels */
  for (l=csidx; l<=ceidx; l++) {
    /* Gather data and error arrays */
    ndat=0; cst=0;
    for (i=0; i<ncomb; i++) {
      if ((lval=l-spec[comb[i][0]].or[comb[i][1]].csidx)>=0 &&
	  (rval=spec[comb[i][0]].or[comb[i][1]].ceidx-l)>=0) {
	if (spec[comb[i][0]].or[comb[i][1]].rdst[lval]==LCLIP ||
	    spec[comb[i][0]].or[comb[i][1]].rdst[lval]==SCLIP)
	  spec[comb[i][0]].or[comb[i][1]].rdst[lval]=1;
	if (spec[comb[i][0]].or[comb[i][1]].rdst[lval]==1) {
	  dat[ndat]=spec[comb[i][0]].or[comb[i][1]].rdfl[lval];
	  err[ndat]=spec[comb[i][0]].or[comb[i][1]].rder[lval];
	  efl[ndat]=spec[comb[i][0]].or[comb[i][1]].rdef[lval];
	  med[ndat]=spec[comb[i][0]].or[comb[i][1]].rdme[lval];
	  res[ndat]=spec[comb[i][0]].or[comb[i][1]].rdres[lval];
	  con[ndat][0]=comb[i][0]; con[ndat][1]=comb[i][1]; con[ndat++][2]=lval;
	} else if (!cst) cst=spec[comb[i][0]].or[comb[i][1]].rdst[lval];
      }
    }

    /* Mark off points which have comparitively large median or statistical
       errors */
    if (ndat) {
      /* Determine maximum value of error and median error for each
	 pixel and then take the minimum value of this over all
	 pixels, multiplied by bigerr, to be the large error
	 threshold */
      bigerr=(MAX(err[0],med[0]));
      for (i=1; i<ndat; i++) bigerr=(MIN(bigerr,(MAX(err[i],med[i]))));
      bigerr*=par->lrgerr;
      i=0; while (i<ndat) {
	if (err[i]>bigerr || med[i]>bigerr) {
	  spec[con[i][0]].or[con[i][1]].rdst[con[i][2]]=LCLIP;
	  j=i; ndat--; while (j<ndat) {
	    dat[j]=dat[j+1]; err[j]=err[j+1]; efl[j]=efl[j+1]; med[j]=med[j+1];
	    res[j]=res[j+1];
	    con[j][0]=con[j+1][0]; con[j][1]=con[j+1][1]; con[j][2]=con[j+1][2];
	    j++;
	  }
	} else i++;
      }
    }

    /* For backwards compatibility with pre-v0.74, check whether
       memory has been freed on alternative arrays or not. We only
       need to do sigma-clipping with error array in next block of
       code (instead of expected fluctuation) when the alternative
       arrays are loaded */
    /* Find the first spectrum with a useful order */
    if (par->version<0.74) {
      for (i=0; i<nspec; i++) {
	for (j=0; j<spec[i].nor; j++) {
	  if (spec[i].or[j].nuse>MINUSE) {
	    if (spec[i].or[j].rdfl_a!=NULL) post074clip=0;
	    break;
    } } } }

    /* Gather final pixel values after sigma-clipping */
    if (!ndat) {
      cspec->fl[l]=0.0; cspec->er[l]=cspec->ef[l]=cspec->res[l]=-INFIN;
      cspec->csq[l]=cspec->ccsq[l]=0.0; cspec->ncb[l]=cspec->nccb[l]=0;
      /* Pre-version 0.66, a distinction was made between cases where
	 no pixels were even available to be combined and where some
	 were available but they were all clipped (for whatever
	 reason). In the first case, the NCLIP flag was applied to the
	 combined spectrum, but in the latter case the status flag of
	 the first available contributing pixel was inherited by the
	 combined spectrum. This was changed in version 0.66 because
	 it was leading to bugs in which the spectra were recombined
	 but were still inheriting the CCLIP flag of the previous
	 combined spectrum. */
      /*
      if (!cst) cspec->st[l]=NCLIP;
      else cspec->st[l]=cst;
      */
      cspec->st[l]=NCLIP;
    } else {
      /* Initialise clip array */
      for (i=0; i<ndat; i++) clip[i]=1;
      /* Gather initial statistics before sigma-clipping */
      if (!stats(dat,err,efl,med,clip,ndat,1,&stat))
	errormsg("UVES_combine_nocont(): Error returned from stats()\n\
\tfor pre-sigma-clip stats for combined pixel %d, wavelength %9.4lf",l+1,
		 cspec->wl[l]);
      cspec->csq[l]=stat.rchisq; cspec->ncb[l]=ndat;
      /* Sigma-clip array (which also does statistics on clipped-filtered
	 array) */
      /* Pre-0.74-version used the statistical uncertainty as the
	 "sigma" in the sigma-clipping algorithm, but this is changed
	 to the expected fluctation in 0.74. This only needs to be
	 done when the alternative arrays the first time and while
	 existing actions are being performed; if their memory has
	 been freed then just use the expected fluctuation from then
	 on. */
      if (!sigclip(dat,err,efl,med,ndat,par->clipsig,&stat,clip,post074clip))
	errormsg("UVES_combine_nocont(): Error doing sigma-clip\n\
\tfor combined spectrum pixel %d, wavelength %9.4lf",l+1,cspec->wl[l]);
      /* Go through clip array and distribute info to source spectra */
      for (i=0,cspec->nccb[l]=0; i<ndat; i++) {
	if (clip[i]==1) (cspec->nccb[l])++;
	else spec[con[i][0]].or[con[i][1]].rdst[con[i][2]]=SCLIP;
      }
      /* Assign final flux and error values to combined spectrum */
      cspec->fl[l]=stat.wmean; cspec->er[l]=stat.ewmean;
      cspec->ef[l]=stat.eflwmean; cspec->ccsq[l]=stat.rchisq; cspec->st[l]=1;
      /* Find the weighted mean resolution of contributing spectra */
      cspec->res[l]=INFIN;
      for (i=0,j=0,sumresonwgtsq=0.0,sum1onwgtsq=0.0; i<ndat; i++) {
	if (clip[i]==1 && res[i]<=0.0) { cspec->res[l]=-INFIN; break; }
	if (clip[i]==1) {
	  wgtsq=med[i]*med[i]; sumresonwgtsq+=res[i]/wgtsq; sum1onwgtsq+=1.0/wgtsq;
	  j++;
	}
      }
      if (cspec->res[l]>0.0) cspec->res[l]=sumresonwgtsq/sum1onwgtsq;
    }
  }

  /* Clean up */
  free(dat); free(err); free(efl); free(med); free(res); free(clip);
  free(*con); free(con);

  return 1;

}
