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

int UVES_rescale_region(spectrum *spec, int nspec, cspectrum *cspec, action *act,
			int nact, double scalclip, double scalerr, params *par) {
  
  double  scale=0.0,escale=0.0,maxsnr=0.0,wgt=0.0,maxwgt=0.0,a=0.0,aerr=0.0,chisq=0.0;
  double  *fit_y=NULL,*fit_e=NULL;
  int     nor=0,noo=0,normax=0,ncov=0,nused=0,region=1,csidx=0,ceidx=0,no=0,nomax=0;
  int     oosidx=0,ooeidx=0,ocsidx=0,oceidx=0;
  int     scsidx=0,sceidx=0;
  int     nfix=0,comb=0,fixed=-1,fit_n=0;
  int     i=0,j=0,k=0,l=0;
  int     *st=NULL;    /* Storage array for original combined status array */
  int     *clip=NULL;  /* Clip array for use with statset */
  int     **fix=NULL;  /* Matrix defining which orders have fixed scaling */
  int     **cov=NULL;  /* Orders covered so far in search of re-scaled orders */ 
  int     **used=NULL; /* Matrix defining which orders have been used so far */ 
  int     **oo=NULL;   /* Matrix defining overalapping orders */
  statset stat;

  /* Find maximum number of orders that could contribute to overlapping
     order matrix */
  for (i=0,nor=0,normax=0; i<nspec; i++) {
    for (j=0,k=0; j<spec[i].nor; j++)
      if (spec[i].or[j].nuse>=MINUSE) { nor++; k++; }
    if (k>normax) normax=k;
  }
  /* Allocate memory for overlapping order, fix, cov and used
     matrices */
  if ((oo=imatrix(nor,2))==NULL)
    errormsg("UVES_rescale_region(): Could not allocate memory\n\
\tfor nor matrix of size %dx%d",nor,2);
  if ((fix=imatrix(nor,2))==NULL)
    errormsg("UVES_rescale_region(): Could not allocate memory\n\
\tfor fix matrix of size %dx%d",nor,2);
  if ((cov=imatrix(nor,2))==NULL)
    errormsg("UVES_rescale_region(): Could not allocate memory\n\
\tfor cov matrix of size %dx%d",nor,2);
  if ((used=imatrix(nspec,normax))==NULL)
    errormsg("UVES_rescale_region(): Could not allocate memory\n\
\tfor used matrix of size %dx%d",nspec,normax);

  /* Go through action array and find which orders have been manually
     scaled by user after most recent recombination */
  for (i=nact-1; i>=0; i--) {
    if (act[i].val && (act[i].act==ARACT || act[i].rcmb==1)) break;
    else if (act[i].val && act[i].act==SOACT) {
      fix[nfix][0]=act[i].i[0]; fix[nfix++][1]=act[i].i[1]; 
    }
  }
  if (!nfix) {
    free(*oo); free(oo); free(*used); free(used); free(*fix); free(fix);
    free(*cov); free(cov); return 2;
  }

  /* Allocate memory and store original combined status array */
  if ((st=iarray(cspec->np))==NULL)
    errormsg("UVES_rescale_region(): Could not allocate memory\n\
\tfor st array of size %d",cspec->np);
  for (i=0; i<cspec->np; i++) st[i]=cspec->st[i];

  /** While there are still uncombined orders, find continous regions and
      combine orders falling in them **/
  nused=0; while (nused<nor) {

    /* Find highest S/N order which hasn't yet been used */
    for (i=0,maxsnr=-INFIN,ncov=0; i<nspec; i++) {
      for (j=0; j<spec[i].nor; j++) {
	if (spec[i].or[j].nuse>=MINUSE && !used[i][j]) {
	  /* Initialize auto-rescale flag for order */
	  spec[i].or[j].ascl=0;
	  if (spec[i].or[j].medsnr>maxsnr) {
	    maxsnr=spec[i].or[j].medsnr; cov[0][0]=i; cov[0][1]=j; ncov=1;
	    csidx=spec[i].or[j].csidx; ceidx=spec[i].or[j].ceidx;
	  }
	}
      }
    } 

    /* If highest SNR order has been fixed then we want to start
       recombining the spectrum, otherwise leave combination flagged
       unchecked */
    for (i=0,scsidx=cspec->np,sceidx=-1,comb=0; i<nfix; i++)
      if (fix[i][0]==cov[0][0] && fix[i][1]==cov[0][1]) {
	scsidx=spec[cov[0][0]].or[cov[0][1]].csidx;
	sceidx=spec[cov[0][0]].or[cov[0][1]].ceidx; comb=1; break;
      }

    region=1; while (region) {
      /* Search for order which most overlaps present combined spectrum */
      for (i=0,nomax=-1,maxwgt=-INFIN; i<nspec; i++) {
	for (j=0; j<spec[i].nor; j++) {
	  if (spec[i].or[j].nuse>=MINUSE) {
	    for (k=0; k<ncov; k++) { if (cov[k][0]==i && cov[k][1]==j) break; }
	    no=(MIN(ceidx,spec[i].or[j].ceidx))-(MAX(csidx,spec[i].or[j].csidx));
	    wgt=spec[i].or[j].medsnr*(double)no;
	    if (k==ncov && no>0 &&
		((!par->scalmeth && no>nomax) || (par->scalmeth==1 && wgt>maxwgt))) {
	      /* If a fixed order was previously found, does the
		 present order overlap with the combination results of
		 that fixed order? */
	      if (comb) {
		if ((MIN(sceidx,spec[i].or[j].ceidx))-
		    (MAX(scsidx,spec[i].or[j].csidx))>=0) {
		  nomax=no; maxwgt=wgt; oo[0][0]=i; oo[0][1]=j;
		}
	      } else { nomax=no; maxwgt=wgt; oo[0][0]=i; oo[0][1]=j; }
	    }
	  }
	}
      }
      /* If no overlapping orders are found then this continous region has
	 been completely combined */
      if (nomax==-1) {
	/* If we have not encountered any orders fixed by the user
	   then copy the coverage array to the used array so that this
	   region is not entered into again */
	if (!comb)
	  for (i=0; i<ncov; i++) { used[cov[i][0]][cov[i][1]]=1; nused++; }
	region=0; break;
      }

      /* Check to see if this order has been fixed by user. Make sure
	 combination flag is set if so. If it wasn't set before then
	 flag this as the first fixed order */
      for (i=0,fixed=-1; i<nfix; i++)
	if (fix[i][0]==oo[0][0] && fix[i][1]==oo[0][1]) {
	  fixed=i; if (!comb) { scsidx=cspec->np; sceidx=-1; }
	  comb=1; break;
	}

      /* If the combination flag is set then carry out the combination */
      if (comb) {
	/* Search for previous orders which overlap with presently
	   selected order */
	for (i=0,l=1; i<nspec; i++) {
	  for (j=0; j<spec[i].nor; j++) {
	    if (spec[i].or[j].nuse>=MINUSE) {
	      for (k=0; k<ncov; k++) {
		if (cov[k][0]==i && cov[k][1]==j &&
		    spec[i].or[j].ceidx>=spec[oo[0][0]].or[oo[0][1]].csidx &&
		    spec[i].or[j].csidx<=spec[oo[0][0]].or[oo[0][1]].ceidx) {
		  oo[l][0]=i; oo[l++][1]=j;
		}
	      }
	    }
	  }
	}
	noo=l;

	/* Find the best-fitting scale factor between this order and
	   the combined spectrum in the overlapping region. There is
	   no need to do this if this order has been fixed by the
	   user */
	if (fixed==-1) {
	  no=nomax;
	  ocsidx=(MAX(csidx,spec[oo[0][0]].or[oo[0][1]].csidx)); oceidx=ocsidx+no;
	  oosidx=ocsidx-spec[oo[0][0]].or[oo[0][1]].csidx; ooeidx=oosidx+no;
	  /* Allocate memory for fitting arrays */
	  if (!par->scalmeth) {
	    if ((fit_y=darray(no))==NULL)
	      errormsg("UVES_rescale_region(): Could not allocate memory\n\
\tfor fit_y array of size %d",no);
	    if ((clip=iarray(no))==NULL)
	      errormsg("UVES_rescale_region(): Could not allocate memory\n\
\tfor clip array of size %d",no);
	  }
	  if ((fit_e=darray(no))==NULL)
	    errormsg("UVES_rescale_region(): Could not allocate memory\n\
\tfor fit_e array of size %d",no);
	  /* As of version 0.40, reinstate pixels in the to-be-scaled
	     order which were previously flagged as having larger
	     errors in a previous combination. */
	  if (par->version>=0.40) {
	    for (i=0,k=oosidx; i<no; i++,k++)
	      if (spec[oo[0][0]].or[oo[0][1]].rdst[k]==LCLIP)
		spec[oo[0][0]].or[oo[0][1]].rdst[k]=1;
	  }
	  /* Choose between scaling methods */
	  scale=1.0; if (!par->scalmeth) {
	    /* Load fitting arrays. Note that a bug has been fixed as
	       of version 0.40 whereby the previous lrgerr status of
	       pixels in the to-be-scaled order is ignored. */
	    for (i=0,j=ocsidx,k=oosidx,fit_n=0; i<no; i++,j++,k++) {
	      if (cspec->st[j]==1 && spec[oo[0][0]].or[oo[0][1]].rdst[k]==1 &&
		  cspec->fl[j]!=0.0) {
		fit_y[fit_n]=spec[oo[0][0]].or[oo[0][1]].rdfl[k]/cspec->fl[j];
		fit_e[fit_n]=sqrt(cspec->er[j]*cspec->er[j]/cspec->fl[j]/
				  cspec->fl[j]+spec[oo[0][0]].or[oo[0][1]].rder[k]*
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
	      /* Do sigma-clipping */
	      if (!sigclip(fit_y,fit_e,NULL,NULL,fit_n,scalclip,&stat,clip))
		errormsg("UVES_rescale_region(): Error returned from sigclip()");
	      scale=stat.wmean; escale=stat.ewmean;
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
			  par->nscalclip,&a,&scale,&aerr,&escale,&chisq,0,1,0))<=0) {
	      if (!i || scale<=0.0 || escale==0.0) scale=1.0;
	    }
	  }
	  free(fit_e);
	  /* Don't scale orders where the scale-factor is poorly
	     determined or negative */
	  if (scale<=0.0 || escale/scale>scalerr) scale=1.0;
	  /* Scale data and errors */
	  for (i=0; i<spec[oo[0][0]].or[oo[0][1]].nrdp; i++) {
	    spec[oo[0][0]].or[oo[0][1]].rdfl[i]/=scale;
	    spec[oo[0][0]].or[oo[0][1]].rder[i]/=scale;
	    spec[oo[0][0]].or[oo[0][1]].rdef[i]/=scale;
	    spec[oo[0][0]].or[oo[0][1]].rdme[i]/=scale;
	  }
	  /* Record scaling for any later undo action */
	  spec[oo[0][0]].or[oo[0][1]].ascl=1;
	  spec[oo[0][0]].or[oo[0][1]].oscl=spec[oo[0][0]].or[oo[0][1]].scl;
	  spec[oo[0][0]].or[oo[0][1]].scl*=scale;
	}

	/* Combine current spectra into combined spectrum */
	if (!UVES_combine_nocont(spec,nspec,cspec,oo,noo,
				 spec[oo[0][0]].or[oo[0][1]].csidx,
				 spec[oo[0][0]].or[oo[0][1]].ceidx,par))
	  errormsg("UVES_rescale_region(): Error attempting to combine\n\
\torders which formed a continous spectral section");

	/* Increase size covered by the orders newly combined after
	   encountering first fixed order */
	scsidx=MIN(scsidx,spec[oo[0][0]].or[oo[0][1]].csidx);
	sceidx=MAX(sceidx,spec[oo[0][0]].or[oo[0][1]].ceidx);

	/* If this is a fixed order, remove it from fixed array */
	if (fixed>=0) {
	  for (i=fixed+1; i<nfix; i++) {
	    fix[i-1][0]=fix[i][0]; fix[i-1][1]=fix[i][1];
	  }
	  nfix--;
	}
      
      }

      /* Register the fact that this order has now been covered */
      cov[ncov][0]=oo[0][0]; cov[ncov++][1]=oo[0][1];

      /* Increase size covered by combined spectrum */
      csidx=MIN(csidx,spec[oo[0][0]].or[oo[0][1]].csidx);
      ceidx=MAX(ceidx,spec[oo[0][0]].or[oo[0][1]].ceidx);
    }
  }

  /* Redefine normalized arrays and return the status array to its
     original state in cases where cosmetic clipping of pixels had
     been done */
  for (i=0; i<cspec->np; i++) {
    if (par->thar<=1) {
      if (st[i]==CCLIP) {
	cspec->st[i]=CCLIP; cspec->no[i]=1.0; cspec->ne[i]=cspec->nf[i]=-INFIN;
      } else {
	cspec->no[i]=cspec->fl[i]/cspec->co[i];
	cspec->ne[i]=cspec->er[i]/cspec->co[i];
	cspec->nf[i]=cspec->ef[i]/cspec->co[i];
      }
    } else {
      if (st[i]==CCLIP) {
	cspec->st[i]=CCLIP; cspec->no[i]=0.0; cspec->ne[i]=cspec->nf[i]=-INFIN;
      } else {
	cspec->no[i]=cspec->fl[i]-cspec->co[i]; cspec->ne[i]=cspec->er[i];
	cspec->ne[i]=cspec->ef[i];
      }
    }
  }

  /* Clean up */
  free(st); free(*oo); free(oo); free(*used); free(used); free(*fix); free(fix);
  free(*cov); free(cov);

  return 1;

}
