/****************************************************************************
* Calculate some basic statistics for the combined spectrum:
* - Wavelength coverage map
* - Median, Min and Max SNR and continuum-to-noise ratio (CNR) at specific wavelengths
* - Resolution in specific wavelengths
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "memory.h"
#include "stats.h"
#include "const.h"
#include "error.h"

int UVES_cspec_stats(spectrum *spec, int nspec, cspectrum *cspec, params *par) {

  double   cdisp=0.0,wls=0.0,wle=0.0;
  double   *snr=NULL,*cnr=NULL,*resnom=NULL,*resarc=NULL;
  double   **wavcov=NULL;
  int      sidx=0,eidx=0,n=0,nval=0;
  int      i=0,j=0,k=0,l=0,m=0;
  statset  stat;
  char     array_name[NAMELEN]="\0";

#define ERR_ARRAY { \
    nferrormsg("UVES_cspec_stats(): Cannot allocate memory for %s array",array_name); \
    return 0; }

  /** Calculate max and median SNR and CNR arrays **/
  /* Determine number of chunks in coarse wavelength grid */
  cdisp=log10(WAVGRIDSIZE/C_C_K+1.0);
  cspec->nc=(int)(log10(WAVGRIDEND/WAVGRIDSTART)/cdisp)+1;
  /* Allocate memory for coarse data arrays */
  sprintf(array_name,"%s","coarse wavel.");
  if ((cspec->cwav=darray(cspec->nc))==NULL) { ERR_ARRAY; }
  sprintf(array_name,"%s","coarse max. SNR");
  if ((cspec->csnrmax=darray(cspec->nc))==NULL) { ERR_ARRAY; }
  sprintf(array_name,"%s","coarse med. SNR");
  if ((cspec->csnrmed=darray(cspec->nc))==NULL) { ERR_ARRAY; }
  if ((cspec->cwav=darray(cspec->nc))==NULL) { ERR_ARRAY; }
  sprintf(array_name,"%s","coarse max. CNR");
  if ((cspec->ccnrmax=darray(cspec->nc))==NULL) { ERR_ARRAY; }
  sprintf(array_name,"%s","coarse med. CNR");
  if ((cspec->ccnrmed=darray(cspec->nc))==NULL) { ERR_ARRAY; }
  sprintf(array_name,"%s","coarse med. nominal resolving power");
  if ((cspec->cresnom=darray(cspec->nc))==NULL) { ERR_ARRAY; }
  sprintf(array_name,"%s","coarse med. arc resolving power");
  if ((cspec->cresarc=darray(cspec->nc))==NULL) { ERR_ARRAY; }

  /* Calculate statistics in each chunk of combined spectrum */
  for (i=0,wle=WAVGRIDSTART,eidx=0; i<cspec->nc; i++) {
    /* Determine starting, central and ending wavelengths for this chunk */
    wls=wle; wle=wls*pow(10.0,cdisp); cspec->cwav[i]=wls*pow(10.0,cdisp*0.5);
    /* Initialize values for this chunk */
    cspec->csnrmax[i]=cspec->csnrmed[i]=cspec->ccnrmax[i]=cspec->csnrmed[i]=0.0;
    cspec->cresnom[i]=cspec->cresarc[i]=0.0;
    /* Determine starting and ending pixels of chunk in combined spectrum */
    if ((sidx=idxdval(&(cspec->wl[eidx]),cspec->np-eidx,wls))==-1) break;
    sidx+=eidx;
    if ((eidx=idxdval(&(cspec->wl[sidx]),cspec->np-sidx,wle))==-1) eidx=cspec->np-1;
    else eidx+=sidx;
    n=eidx-sidx+1;
    /* Allocate memory for data arrays */
    sprintf(array_name,"%s","snr"); if ((snr=darray(n))==NULL) { ERR_ARRAY; }
    sprintf(array_name,"%s","cnr"); if ((cnr=darray(n))==NULL) { ERR_ARRAY; }
    sprintf(array_name,"%s","resnom"); if ((resnom=darray(nspec))==NULL) { ERR_ARRAY; }
    sprintf(array_name,"%s","resarc"); if ((resarc=darray(nspec))==NULL) { ERR_ARRAY; }
    /* Calculate SNR and CNR arrays */
    for (j=sidx,k=0; j<=eidx; j++) {
      if (cspec->st[j]==1 && cspec->wl[j]>=wls && cspec->wl[j]<=wle) {
	snr[k]=cspec->no[j]/cspec->ne[j]; cnr[k++]=1.0/cspec->ne[j];
      }
    }
    nval=k;
    if (nval) {
      /* Determine max. SNR and CNR */
      cspec->csnrmax[i]=djmax(snr,nval); cspec->ccnrmax[i]=djmax(cnr,nval);
      /* Determine med. SNR and CNR */
      stat.med=0.0; if (!median(snr,NULL,nval,&stat,0))
	nferrormsg("UVES_cspec_stats(): Error calculating median SNR for\n\
\tregion %lf-%lf of combined spectrum which contains %d valid pixels",wls,wle,nval);
      cspec->csnrmed[i]=stat.med;
      stat.med=0.0; if (!median(cnr,NULL,nval,&stat,0))
	nferrormsg("UVES_cspec_stats(): Error calculating median CNR for\n\
\tregion %lf-%lf of combined spectrum which contains %d valid pixels",wls,wle,nval);
      cspec->ccnrmed[i]=stat.med;
      /* Calculate the median resolving power using the hard-coded,
	 slit-width--resolution product values for each arm and binning
	 in UVES_popler.h */
      /* Allocate memory for data array */
      /* Fill res arrays with nominal and arc resolving powers of
	 spectra based on their slit-widths and actual arc extractions */
      for (j=0,k=0,n=0; j<nspec; j++) {
	/* Test whether this spectrum contributes to this chunk simply
	   by seeing whether the centre of this chunk lies within
	   total wavelength range of spectrum. This ignores whether
	   the orders were included in the final combination of the
	   spectrum, but it's simple */
	l=0; while (l<spec[j].nor && spec[j].or[l].nuse<MINUSE) l++;
	m=spec[j].nor-1; while (m>=0 && spec[j].or[m].nuse<MINUSE) m--;
	if (l<spec[j].nor && m>=0 && l<=m &&
	    cspec->cwav[i]>=cspec->wl[spec[j].or[l].csidx] &&
	    cspec->cwav[i]<=cspec->wl[spec[j].or[m].ceidx]) {
	  l=(spec[j].binx==1) ? 0 : 1;
	  m=(spec[j].cwl<UVESBORR) ? 0 : 1;
	  if (spec[j].sw>0.0)
	    resnom[k++]=UVES_NOM_RES[l][m][0]/spec[j].sw+UVES_NOM_RES[l][m][1]+
	      UVES_NOM_RES[l][m][2]*spec[j].sw;
	  if (spec[j].arcfwhm>0.0) resarc[n++]=C_C_K/spec[j].arcfwhm;
	}
      }
      if (k>0) {
	stat.med=0.0; if (!median(resnom,NULL,k,&stat,0))
	  nferrormsg("UVES_cspec_stats(): Error calculating median nominal\n\
\tresolving power for region %lf-%lf of combined spectrum",wls,wle);
	cspec->cresnom[i]=stat.med;
      }
      if (n>0) {
	stat.med=0.0; if (!median(resarc,NULL,n,&stat,0))
	  nferrormsg("UVES_cspec_stats(): Error calculating median arc\n\
\tresolving power for region %lf-%lf of combined spectrum",wls,wle);
	cspec->cresarc[i]=stat.med;
      }
    }
    /* Free memory for data arrays */
    free(snr); free(cnr); free(resnom); free(resarc);
  }

  /** Determine wavelength coverage map **/
  /* Allocate memory for temporary wavelength coverage map */
  if ((wavcov=dmatrix(cspec->np/2+1,4))==NULL) {
    nferrormsg("UVES_cspec_stats(): Cannot allocate memory for temporary\n\
\twavelength coverage matrix of size %dx%d",cspec->np/2+1,4); \
    return 0;
  }
  /* Find pairs of regions of valid and invalid pixels. Record the
     starting and ending wavelengths of the valid regions and the
     sizes (in km/s) of the valid and invalid regions */
  i=n=0;
  while (i<cspec->np) {
    while (i<cspec->np && cspec->st[i]!=1) i++; wavcov[n][0]=cspec->wl[i];
    j=i; while (j<cspec->np && cspec->st[j]==1) j++; j--; wavcov[n][1]=cspec->wl[j];
    wavcov[n][2]=2.0*C_C_K*(wavcov[n][1]-wavcov[n][0])/(wavcov[n][1]+wavcov[n][0]);
    k=j+1; while (k<cspec->np && cspec->st[k]!=1) k++;
    wavcov[n++][3]=2.0*C_C_K*(cspec->wl[k-1]-cspec->wl[j])/(cspec->wl[k-1]+cspec->wl[j]);
    i=k;
  }
  /* Go through temporary coverage map and decide which regions are
     large enough and which gaps are small enough. Delete any regions
     that aren't large enough, or are separated from others (that are
     large enough) by small enough gaps, and shuffle the remaining
     regions up the matrix */
  cspec->nwavcov=0;
  i=0; while (i<n) {
    /* Find next gap that is too large */
    j=i; while (j<n && wavcov[j][3]<=WAVCOVMAXGAP) j++;
    /* If needed, edit the coverage matrix */
    if (j>i && j<n) {
      wavcov[i][1]=wavcov[j][1];
      wavcov[i][2]=2.0*C_C_K*(wavcov[i][1]-wavcov[i][0])/(wavcov[i][1]+wavcov[i][0]);
      wavcov[i][3]=wavcov[j][3];
      for (k=i+1; k<n; k++) { for (l=0; l<4; l++) wavcov[k][l]=wavcov[j+k-i][l]; }
      n-=j-i;
    }
    /* If this region is big enough then count it */
    if (wavcov[i][2]>WAVCOVMINREG) cspec->nwavcov++;
    i++;
  }
  /* Allocate memory for permanent wavelength coverage map */
  if ((cspec->wavcov=dmatrix(cspec->nwavcov,2))==NULL) {
    nferrormsg("UVES_cspec_stats(): Cannot allocate memory for wavelength\n\
\tcoverage matrix of size %dx%d",cspec->nwavcov,2); \
    return 0;
  }
  /* Copy large enough wavelength coverage regions to permanent map */
  for (i=0,j=0; i<n; i++) {
    if (wavcov[i][2]>WAVCOVMINREG) {
      cspec->wavcov[j][0]=wavcov[i][0]; cspec->wavcov[j][1]=wavcov[i][1];
      j++;
    }
  }
  /* Free temporary matrix */
  free(*wavcov); free(wavcov);

  /* Find the total number of exposures and total exposure time. At
     the moment, this is treated straight-forwardly, except for UVES
     because it has two arms with one chip in the blue arm and 3 chips
     in the red arm. Depending on whether a dichroic exposure was
     taken, or not, there can be 1, 2 or 3 files per exposure. At the
     moment, it is just assumed that there was one exposure per
     observation block (OB) for UVES so that the "HIERARCH ESO OBS ID"
     header card can be used to identify files originating from the
     same exposure. */
  for (i=0,cspec->nexp=0,cspec->texp=0.0; i<nspec; i++) {
    if (spec[i].ftype!=FTUVES) { cspec->nexp++; cspec->texp+=spec[i].etime; }
    else {
      for (j=i+1; j<nspec; j++) {
	if (spec[j].ftype==FTUVES && spec[j].obid==spec[i].obid &&
	    fabs(spec[j].jd-spec[i].jd)*C_SECDAY<0.1*(MIN(spec[i].etime,spec[j].etime))) {
	  if (spec[i].etime>spec[j].etime) cspec->texp+=spec[i].etime-spec[j].etime;
	  break;
	}
      }
      if (j==nspec) { cspec->nexp++; cspec->texp+=spec[i].etime; }
    }
  }

  return 1;

}
