/****************************************************************************
* Combine all orders of all exposures
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "stats.h"
#include "memory.h"
#include "error.h"

int UVES_combine_spec(spectrum *spec, int nspec, cspectrum *cspec,
		      params *par) {
  
  double   bigerr=0.0;
  double   *dat=NULL,*err=NULL,*efl=NULL,*med=NULL;
  int      ndat=0,maxndat=0,cst=0;
  int      lval=0,rval=0;
  int      i=0,j=0,l=0;
  int      *clip=NULL,**con=NULL;
  statset  stat;

  /* Find maximum number of possible contributions to continuum at each
     point in combined spectrum and allocate memory to data arrays */
  for (i=0; i<nspec; i++) {
    for (j=0; j<spec[i].nor; j++) maxndat++; 
  }
  if ((dat=darray(maxndat))==NULL)
    errormsg("UVES_combine_spec(): Could not allocate\n\
\tmemory for dat array of size %d",maxndat);
  if ((err=darray(maxndat))==NULL)
    errormsg("UVES_combine_spec(): Could not allocate\n\
\tmemory for err array of size %d",maxndat);
  if ((efl=darray(maxndat))==NULL)
    errormsg("UVES_combine_spec(): Could not allocate\n\
\tmemory for efl array of size %d",maxndat);
  if ((med=darray(maxndat))==NULL)
    errormsg("UVES_combine_spec(): Could not allocate\n\
\tmemory for med array of size %d",maxndat);
  if ((clip=iarray(maxndat))==NULL)
    errormsg("UVES_combine_spec(): Could not allocate\n\
\tmemory for clip array of size %d",maxndat);
  if ((con=imatrix(maxndat,3))==NULL)
    errormsg("UVES_combine_spec(): Could not allocate\n\
\tmemory for con matrix of size %dx%d",maxndat,3);

  for (l=0; l<cspec->np; l++) {
    /* Gather data and error arrays */
    ndat=0; cst=0;
    for (i=0; i<nspec; i++) {
      if (spec[i].comb) {
	for (j=0; j<spec[i].nor; j++) {
	  if (spec[i].or[j].nuse>MINUSE) {
	    if ((lval=l-spec[i].or[j].csidx)>=0 &&
		(rval=spec[i].or[j].ceidx-l)>=0) {
	      if (spec[i].or[j].rdst[lval]==LCLIP ||
		  spec[i].or[j].rdst[lval]==SCLIP) spec[i].or[j].rdst[lval]=1;
	      if (spec[i].or[j].rdst[lval]==1) {
		dat[ndat]=spec[i].or[j].rdfl[lval];
		err[ndat]=spec[i].or[j].rder[lval];
		efl[ndat]=spec[i].or[j].rdef[lval];
		med[ndat]=spec[i].or[j].rdme[lval];
		con[ndat][0]=i; con[ndat][1]=j; con[ndat++][2]=lval;
	      }
	      else if (!cst) cst=spec[i].or[j].rdst[lval];
	    }
	  }
	}
      }
    }
    
    /* Mark off points which have comparitively large errors */
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
	    con[j][0]=con[j+1][0]; con[j][1]=con[j+1][1]; con[j][2]=con[j+1][2];
	    j++;
	  }
	}
	else i++;
      }
    }

    /* Gather final pixel values after sigma-clipping */
    if (!ndat) {
      cspec->fl[l]=cspec->co[l]; cspec->no[l]=1.0;
      cspec->er[l]=cspec->ef[l]=-INFIN; cspec->ne[l]=cspec->nf[l]=-INFIN;
      cspec->csq[l]=cspec->ccsq[l]=0.0; cspec->ncb[l]=cspec->nccb[l]=0;
      if (!cst) cspec->st[l]=NCLIP;
      else cspec->st[l]=cst;
    }
    else {
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
      if (!sigclip(dat,err,efl,med,ndat,par->clipsig,&stat,clip))
	errormsg("UVES_combine_nocont(): Error doing sigma-clip\n\
\tfor combined spectrum pixel %d, wavelength %9.4lf",l+1,cspec->wl[l]);
      /* Go through clip array and distribute info to source spectra */
      for (i=0,cspec->nccb[l]=0; i<ndat; i++) {
	if (clip[i]==1) (cspec->nccb[l])++;
	else spec[con[i][0]].or[con[i][1]].rdst[con[i][2]]=SCLIP;
      }
      /* Assign final flux and error values to combined spectrum */
      cspec->fl[l]=stat.wmean; cspec->er[l]=stat.ewmean;
      cspec->ef[l]=stat.eflwmean; cspec->ccsq[l]=stat.rchisq;
      if (par->thar<=1) {
	if (cspec->st[l]==CCLIP) {	
	  cspec->no[l]=1.0; cspec->ne[l]=cspec->nf[l]=-INFIN;
	} else {
	  cspec->no[l]=cspec->fl[l]/cspec->co[l];
	  cspec->ne[l]=cspec->er[l]/cspec->co[l];
	  cspec->nf[l]=cspec->ef[l]/cspec->co[l]; cspec->st[l]=1;
	}
      } else {
	if (cspec->st[l]==CCLIP) {
	  cspec->no[l]=0.0; cspec->ne[l]=cspec->nf[l]=-INFIN;
	} else {
	  cspec->no[l]=cspec->fl[l]-cspec->co[l]; cspec->ne[l]=cspec->er[l];
	  cspec->nf[l]=cspec->ef[l]; cspec->st[l]=1;
	}
      }
    }
  }

  /* Since UVES_set_wavelen_scale adds an extra 1 or 2 pixels on the end,
     just in case there is some residual flux there, here we should check
     these pixels' flux and reject them if there's none */
  if (cspec->fl[cspec->np-1]<DRNDTOL) {
    cspec->np--; if (cspec->fl[cspec->np-1]<DRNDTOL) cspec->np--;
  }

  /* Clean up */
  free(dat); free(err); free(efl); free(clip); free(*con); free(con);

  return 1;

}
