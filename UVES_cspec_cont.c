/****************************************************************************
* Find a continuum for each contigious region of combined spectrum and
* apply that to the orders underlying it.
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "stats.h"
#include "const.h"
#include "error.h"

int UVES_cspec_cont(spectrum *spec, int nspec, cspectrum *cspec, params *par) {

  double   grad=0.0;
  int      sidx=0,eidx=0,lyaidx=0,beidx=0,rsidx=0,ssecidx=0,esecidx=0;
  int      ncont=0;
  int      i=0,j=0,k=0;
  statset  stat;

  /* Do not fit a continuum if we are dealing with ThAr spectra or if
     we've been asked not to by user */
  if (par->thar<=1 && !par->nocont) {
    /* Find first and last valid pixel in cspec */
    i=0; while (cspec->st[i]<1) i++; sidx=i;
    i=cspec->np-1; while (cspec->st[i]<1) i--; eidx=i;
    /* Find the index of (vc below) Lyman alpha line */
    if (fabs(par->zem)>DRNDTOL) 
      lyaidx=idxdval(cspec->wl,cspec->np,LYA*(1.0-par->vlya/C_C_K)*(1.0+par->zem));
    else lyaidx=-2;
    if (lyaidx>=0) lyaidx=(MAX(lyaidx,sidx));
  
    /*** Fit contigious regions below LyA emission line of QSO ***/
    if (lyaidx==-1 || lyaidx>sidx+MINUSE) {
      /* Find last valid pixel before Lya marker */
      beidx=lyaidx-1; while (beidx>sidx+MINUSE && cspec->st[beidx]<1) beidx--;
      /* Begin loop to fit contigious regions of continuum */
      ssecidx=sidx; esecidx=ssecidx+1;
      while (esecidx<beidx) {
	/* Determine start and end of contigious region */
	while (esecidx<=beidx && cspec->st[esecidx]!=RCLIP &&
	       cspec->st[esecidx]!=NCLIP) esecidx++; esecidx--;
	ncont=esecidx-ssecidx+1;
	/* Fit continuum over contigious region */
	if (ncont>MINUSE) {
	  if (!UVES_chunk_cont(&(cspec->wl[ssecidx]),&(cspec->fl[ssecidx]),
			       &(cspec->er[ssecidx]),&(cspec->co[ssecidx]),
			       &(cspec->st[ssecidx]),ncont,par->vclya,
			       par->rsiglyal,par->rsiglyau,par->pctllya,
			       par->ftyplya,par->cordlya))
	    errormsg("UVES_cspec_cont(): Error fitting continuum to section\n\
\tof combined spectrum from wavelength %lf to %lf",
		     cspec->wl[ssecidx],cspec->wl[esecidx]);
	}
	else {
	  if (!median(&(cspec->fl[ssecidx]),NULL,ncont,&stat,0))
	    errormsg("UVES_cspec_cont(): Error returned from median() when\n\
\ttaking median of combined spectrum from wavelength %lf to %lf",
		     cspec->wl[ssecidx],cspec->wl[esecidx]);
	  for (i=ssecidx; i<=esecidx; i++) cspec->co[i]=stat.med;
	}
	/* Find start of next region and fill in continuum for intervening
	   pixels */
	ssecidx=esecidx+1;
	while (ssecidx<beidx && (cspec->st[ssecidx]==RCLIP ||
				 cspec->st[ssecidx]==NCLIP)) ssecidx++;
	grad=(cspec->co[ssecidx]-cspec->co[esecidx])/(double)(ssecidx-esecidx);
	for (i=esecidx+1; i<ssecidx; i++)
	  cspec->co[i]=cspec->co[esecidx]+grad*(double)(i-esecidx);
	esecidx=ssecidx+1;
      }
    }

    /*** Fit contigious regions above LyA emission line of QSO ***/
    if (lyaidx==-2 || eidx-lyaidx>MINUSE) {
      /* Find first valid pixel after Lya marker */
      rsidx=(MAX(lyaidx,0));
      while (eidx-rsidx>MINUSE && cspec->st[rsidx]<1) rsidx++;
      /* Begin loop to fit contigious regions of continuum */
      ssecidx=rsidx; esecidx=ssecidx+1;
      while (esecidx<eidx) {
	/* Determine start and end of contigious region */
	while (esecidx<=eidx && cspec->st[esecidx]!=RCLIP &&
	       cspec->st[esecidx]!=NCLIP) esecidx++; esecidx--;
	ncont=esecidx-ssecidx+1;
	/* Fit continuum over contigious region */
	if (ncont>MINUSE) {
	  if (!UVES_chunk_cont(&(cspec->wl[ssecidx]),&(cspec->fl[ssecidx]),
			       &(cspec->er[ssecidx]),&(cspec->co[ssecidx]),
			       &(cspec->st[ssecidx]),ncont,par->vcred,
			       par->rsigredl,par->rsigredu,par->pctlred,
			       par->ftypred,par->cordred))
	    errormsg("UVES_cspec_cont(): Error fitting continuum to section\n\
\tof combined spectrum from wavelength %lf to %lf",
		     cspec->wl[ssecidx],cspec->wl[esecidx]);
	}
	else {
	  if (!median(&(cspec->fl[ssecidx]),NULL,ncont,&stat,0))
	    errormsg("UVES_cspec_cont(): Error returned from median() when\n\
\ttaking median of combined spectrum from wavelength %lf to %lf",
		     cspec->wl[ssecidx],cspec->wl[esecidx]);
	  for (i=ssecidx; i<=esecidx; i++) cspec->co[i]=stat.med;
	}
	/* Find start of next region and fill in continuum for intervening
	   pixels */
	ssecidx=esecidx+1;
	while (ssecidx<eidx && (cspec->st[ssecidx]==RCLIP ||
				cspec->st[ssecidx]==NCLIP)) ssecidx++;
	esecidx=ssecidx+1;
      }
    }

    /* Go through combined spectrum and fill in gaps */
    sidx=eidx=0; while (sidx<cspec->np) {
      while (sidx<cspec->np && cspec->st[sidx]!=RCLIP && cspec->st[sidx]!=NCLIP)
	sidx++; eidx=sidx;
      while (eidx<cspec->np && (cspec->st[eidx]==RCLIP || cspec->st[eidx]==NCLIP))
	eidx++; eidx--;
      if (sidx==cspec->np) break;
      else if (!sidx) {
	for (i=sidx; i<=eidx; i++) {
	  cspec->co[i]=cspec->fl[i]=cspec->co[eidx+1];
	  cspec->er[i]=cspec->ef[i]=-INFIN;
	}
      }
      else if (eidx==cspec->np-1) {
	for (i=sidx; i<=eidx; i++) {
	  cspec->co[i]=cspec->fl[i]=cspec->co[sidx-1];
	  cspec->er[i]=cspec->ef[i]=-INFIN;
	}
      }
      else {
	grad=(cspec->co[eidx+1]-cspec->co[sidx-1])/(double)(eidx-sidx+2);
	for (i=sidx; i<=eidx; i++) {
	  cspec->co[i]=cspec->fl[i]=cspec->co[sidx-1]+grad*(double)(i+1-sidx);
	  cspec->er[i]=cspec->ef[i]=-INFIN;
	}
      }
      sidx=eidx+1;
    }
  }

  /* Set the normalized flux and error arrays */
  if (par->thar<=1) {
    if (par->nocont==CNTNON) for (i=0; i<cspec->np; i++) cspec->co[i]=1.0;
    for (i=0; i<cspec->np; i++) {
      if (cspec->st[i]==1) {
	cspec->no[i]=cspec->fl[i]/cspec->co[i];
	cspec->ne[i]=cspec->er[i]/cspec->co[i];
	cspec->nf[i]=cspec->ef[i]/cspec->co[i];
      } else {
	cspec->no[i]=1.0; cspec->ne[i]=cspec->nf[i]=-INFIN;
      }
    }
  }
  else {
    for (i=0; i<cspec->np; i++) {
      cspec->co[i]=0.0; cspec->no[i]=cspec->fl[i]; cspec->ne[i]=cspec->er[i];
      cspec->nf[i]=cspec->ef[i];
    }
  }

  /* Loop over spectra to set order continuua equal to combined continua */
  for (i=0; i<nspec; i++) {
    /* Loop over useful orders ... */
    for (j=0; j<spec[i].nor; j++) {
      if (spec[i].or[j].nuse>=MINUSE) {
	for (k=0; k<spec[i].or[j].nrdp; k++)
	  spec[i].or[j].rdco[k]=cspec->co[spec[i].or[j].csidx+k];
      }
    }
  }
  
  return 1;

}
