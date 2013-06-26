/****************************************************************************
* Combined the spectra by finding continuous regions of spectrum and
* hierachically combining nearby spectra, matching only the overlapping
* regions in the flux domain.
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "stats.h"
#include "fit.h"
#include "memory.h"
#include "error.h"

int UVES_combine_region(spectrum *spec, int nspec, cspectrum *cspec, params *par) {

  double  scale=0.0,scalerr=0.0,maxsnr=0.0,wgt=0.0,maxwgt=0.0,a=0.0,aerr=0.0,chisq=0.0;
  double  *fit_y=NULL,*fit_e=NULL;
  int     nor=0,noo=0,normax=0,nused=0,region=1,csidx=0,ceidx=0,no=0,nomax=0;
  int     oosidx=0,ooeidx=0,ocsidx=0,oceidx=0;
  int     rank=0;
  int     fit_n=0;
  int     i=0,j=0,k=0,l=0;
  int     *clip=NULL;  /* Clip array for use with statset */
  int     **used=NULL; /* Matrix defining which orders have been used so far */ 
  int     **oo=NULL;   /* Matrix defining overalapping orders */
  statset stat;

  /* Scale each spectrum+continuum according to its exposure time and 
     Find maximum number of orders that could contribute to overlapping
     order matrix */
  for (i=0,nor=0,normax=0; i<nspec; i++) {
    scale=(par->thar==2 || par->rankspec) ? 1.0 : 3600.0/spec[i].etime;
    // scale=(par->rankspec) ? 1.0 : 3600.0/spec[i].etime;
    for (j=0,k=0; j<spec[i].nor; j++) {
      /* Check first to see if order is useful */
      if (spec[i].or[j].nuse>=MINUSE) {
	nor++; k++;
	if (par->thar<=1) {
	  spec[i].or[j].scl=scale;
	  for (l=0; l<spec[i].or[j].nrdp; l++) {
	    spec[i].or[j].rdfl[l]*=scale; spec[i].or[j].rdme[l]*=scale;
	    if (spec[i].or[j].rdst[l]==1) {
	      spec[i].or[j].rder[l]*=scale; spec[i].or[j].rdef[l]*=scale;
	    }
	  }
	}
      }
    }
    if (k>normax) normax=k;
  }
  rank=nor;

  /* Allocate memory for overlapping order matrix and "used" matrix */
  if ((oo=imatrix(nor,2))==NULL)
    errormsg("UVES_combine_region(): Could not allocate memory\n\
\tfor nor matrix of size %dx%d",nor,2);
  if ((used=imatrix(nspec,normax))==NULL)
    errormsg("UVES_combine_region(): Could not allocate memory\n\
\tfor used matrix of size %dx%d",nspec,normax);

  /** While there are still uncombined orders, find continous regions and
      combine orders falling in them **/
  nused=0; while (nused<nor) {
    /* Find highest S/N order which hasn't yet been used */
    if (!par->combmeth && par->rankspec && !nused) {
      for (j=0,maxsnr=-INFIN; j<spec[0].nor; j++) {
	if (spec[0].or[j].nuse>=MINUSE && spec[0].or[j].medsnr>maxsnr) {
	  maxsnr=spec[0].or[j].medsnr; oo[0][0]=0; oo[0][1]=j;
	}
      }
    } else {
      for (i=0,maxsnr=-INFIN; i<nspec; i++) {
	for (j=0; j<spec[i].nor; j++) {
	  if (spec[i].or[j].nuse>=MINUSE && !used[i][j] &&
	      spec[i].or[j].medsnr>maxsnr) {
	    maxsnr=spec[i].or[j].medsnr; oo[0][0]=i; oo[0][1]=j;
	  }
	}
      }
    }
    /** Copy this spectrum to the combined spectrum **/
    csidx=spec[oo[0][0]].or[oo[0][1]].csidx; ceidx=spec[oo[0][0]].or[oo[0][1]].ceidx;
    for (i=0,j=csidx; i<spec[oo[0][0]].or[oo[0][1]].nrdp; i++,j++) {
      cspec->fl[j]=spec[oo[0][0]].or[oo[0][1]].rdfl[i];
      cspec->er[j]=spec[oo[0][0]].or[oo[0][1]].rder[i];
      cspec->ef[j]=spec[oo[0][0]].or[oo[0][1]].rdef[i];
      cspec->res[j]=spec[oo[0][0]].or[oo[0][1]].rdres[i];
      cspec->st[j]=spec[oo[0][0]].or[oo[0][1]].rdst[i]; cspec->ncb[j]=cspec->nccb[j]=1;
    }
    /* Record the use of the high-snr order and rank it */
    used[oo[0][0]][oo[0][1]]=1; nused++;
    spec[oo[0][0]].or[oo[0][1]].crank=rank--;
    
    /** Algorithm: combine the high-snr order with the order which most
        overlaps with it **/
    region=1; while (region) {
      /* Search for order which most overlaps present combined spectrum */
      for (i=0,nomax=-1,maxwgt=-INFIN; i<nspec; i++) {
	for (j=0; j<spec[i].nor; j++) {
	  no=(MIN(ceidx,spec[i].or[j].ceidx))-(MAX(csidx,spec[i].or[j].csidx));
	  wgt=spec[i].or[j].medsnr*(double)no;
	  if (spec[i].or[j].nuse>=MINUSE && !used[i][j] && no>0 &&
	      ((!par->scalmeth && no>nomax) || (par->scalmeth==1 && wgt>maxwgt))) {
	    nomax=no; maxwgt=wgt; oo[0][0]=i; oo[0][1]=j; 
	  }
	}
      }
      /* If no overlapping orders are found then this continous region has
	 been completely combined */
      if (nomax==-1) { region=0; break; }

      /* Search for already-combined orders which overlap current newly
	 selected order */
      for (i=0,l=1; i<nspec; i++) {
	for (j=0; j<spec[i].nor; j++) {
	  if (spec[i].or[j].nuse>=MINUSE && used[i][j] &&
	      spec[i].or[j].ceidx>=spec[oo[0][0]].or[oo[0][1]].csidx &&
	      spec[i].or[j].csidx<=spec[oo[0][0]].or[oo[0][1]].ceidx) {
	    oo[l][0]=i; oo[l++][1]=j;
	  }
	}
      }
      noo=l;

      /*** Find scale factor relating combined spectrum and newly found
	   overlapping order ***/
      no=nomax;
      ocsidx=(MAX(csidx,spec[oo[0][0]].or[oo[0][1]].csidx)); oceidx=ocsidx+no;
      oosidx=ocsidx-spec[oo[0][0]].or[oo[0][1]].csidx; ooeidx=oosidx+no;
      if (par->thar<=1) {
	/** For object frames, just find best fitting ratio of overalapping
	    orders **/
	/* Allocate memory for fitting arrays */
	if (!par->scalmeth) {
	  if ((fit_y=darray(no))==NULL)
	    errormsg("UVES_combine_region(): Could not allocate memory\n\
\tfor fit_y array of size %d",no);
	  if ((clip=iarray(no))==NULL)
	    errormsg("UVES_combine_region(): Could not allocate memory\n\
\tfor clip array of size %d",no);
	}
	if ((fit_e=darray(no))==NULL)
	  errormsg("UVES_combine_region(): Could not allocate memory\n\
\tfor fit_e array of size %d",no);
	/* Choose between scaling methods */
	scale=1.0; if (!par->scalmeth) {
	  /* Load fitting arrays */
	  for (i=0,j=ocsidx,k=oosidx,fit_n=0; i<no; i++,j++,k++) {
	    if (cspec->st[j]==1 && spec[oo[0][0]].or[oo[0][1]].rdst[k]==1 &&
		cspec->fl[j]!=0.0) {
	      fit_y[fit_n]=spec[oo[0][0]].or[oo[0][1]].rdfl[k]/cspec->fl[j];
	      fit_e[fit_n]=sqrt(cspec->er[j]*cspec->er[j]/cspec->fl[j]/cspec->fl[j]+
				spec[oo[0][0]].or[oo[0][1]].rder[k]*
				spec[oo[0][0]].or[oo[0][1]].rder[k]/
				spec[oo[0][0]].or[oo[0][1]].rdfl[k]/
				spec[oo[0][0]].or[oo[0][1]].rdfl[k]);
	      fit_n++;
	    }
	  }
	  /* Sigma-clip array iteratively */
	  if (fit_n>1) {
	    /* Initialize clip array */
	    for (i=0; i<fit_n; i++) clip[i]=1;
	    /* Do the sigma-clipping */
	    if (!sigclip(fit_y,fit_e,NULL,NULL,fit_n,par->scalclip,&stat,clip))
	      errormsg("UVES_combine_region(): Error returned from sigclip()");
	    scale=stat.wmean; scalerr=stat.ewmean;
	  } else scale=1.0;
	  free(fit_y); free(clip);
	} else {
	  /* Load maximum error array for fitting */
	  for (i=0,j=oosidx; i<no; i++,j++)
	    fit_e[i]=(MAX(spec[oo[0][0]].or[oo[0][1]].rder[j],
			  spec[oo[0][0]].or[oo[0][1]].rdme[j]));
	  /* Determine best-fit ratio */
	  if ((i=fitexy(&(cspec->fl[ocsidx]),
			 &(spec[oo[0][0]].or[oo[0][1]].rdfl[oosidx]),
			 &(cspec->er[ocsidx]),fit_e,&(cspec->st[ocsidx]),
			 &(spec[oo[0][0]].or[oo[0][1]].rdst[oosidx]),no,par->scalclip,
			par->nscalclip,&a,&scale,&aerr,&scalerr,&chisq,0,1,0))<=0) {
	    if (!i || scale<=0.0 || scalerr==0.0) scale=1.0;
	  }
	}
	free(fit_e);
	/* Don't scale orders where the scale-factor is poorly
	   determined or negative */
	if (scale<=0.0 || scalerr/scale>par->scalerr) scale=1.0;
	/* Scale data and errors */
	for (i=0; i<spec[oo[0][0]].or[oo[0][1]].nrdp; i++) {
	  spec[oo[0][0]].or[oo[0][1]].rdfl[i]/=scale;
	  spec[oo[0][0]].or[oo[0][1]].rder[i]/=scale;
	  spec[oo[0][0]].or[oo[0][1]].rdef[i]/=scale;
	  spec[oo[0][0]].or[oo[0][1]].rdme[i]/=scale;
	}
	/* Update scale factor applied to this order */
	spec[oo[0][0]].or[oo[0][1]].scl/=scale;
      }

      /** Combine current spectra (making up continuous region) into a
	  combined spectrum **/
      /* Do the combination */
      if (!UVES_combine_nocont(spec,nspec,cspec,oo,noo,
	  spec[oo[0][0]].or[oo[0][1]].csidx,spec[oo[0][0]].or[oo[0][1]].ceidx,par))
	errormsg("UVES_combine_region(): Error attempting to combine\n\
\torders which formed a continous spectral section");
      /* Record use of this order, rank it and move on to next
	 available overlapping section */
      used[oo[0][0]][oo[0][1]]=1; nused++; spec[oo[0][0]].or[oo[0][1]].crank=rank--;

      /* Define new starting and ending pixels */
      csidx=(MIN(csidx,spec[oo[0][0]].or[oo[0][1]].csidx));
      ceidx=(MAX(ceidx,spec[oo[0][0]].or[oo[0][1]].ceidx));
    }
  }

  /* Clean up */
  free(*used); free(used); free(*oo); free(oo);

  return 1;

}
