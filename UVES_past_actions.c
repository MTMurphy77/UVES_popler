/****************************************************************************
* Execute all the historical actions that have been stored in the action
* array and which have been read in from an old UPL file
*
* opt = 0  : Execute all actions in array
* opt = n  : Execute only first n actions in array
* opt = -n : Execute only last n actions in array
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "fit.h"
#include "memory.h"
#include "error.h"

int UVES_past_actions(spectrum *spec, int nspec, cspectrum *cspec, action *act,
		      int nact, int opt, params *par) {

  double    frac=0.0,pctl=0.0,scale=0.0;
  double    *fco=NULL,*fx=NULL,*fy=NULL,*fy2=NULL;
  int       op1=0,op2=0,op3=0,op4=0;
  int       cp1=0,cp2=0,cp3=0,cp4=0;
  int       fp1=0,fp2=0,fnp=0,frs=0,fbsnp=0,fbenp=0;
  int       sact=0,eact=0,ncon=0,idx=0,status=0;
  int       i=0,j=0,k=0,l=0;
  int       *fst=NULL;
  int       **con=NULL;

  /* Allocate memory for contribution matrix */
  for (i=0,k=0; i<nspec; i++) {
    for (j=0; j<spec[i].nor; j++) {
      if (spec[i].or[j].nuse>=MINUSE) k++;
    }
  }
  if ((con=imatrix(k,4))==NULL)
    errormsg("UVES_past_actions(): Cannot allocate memory\n\
\tfor con matrix of size %dx%d",k,4);

  /* Determine which action to start on and how many to execute */
  if (!opt) { sact=0; eact=nact-1; }
  else if (opt>0) {
    if (opt>nact) {
      nferrormsg("UVES_past_actions(): Number of actions to\n\
\texecute (opt=%d) is greater than number available (=%d)",opt,nact); return 0;
    }
    sact=0; eact=opt-1;
  }
  else {
    if (-1*opt>nact) {
      nferrormsg("UVES_past_actions(): Number of actions to\n\
\texecute (opt=%d) is greater than number available (=%d)",-1*opt,nact); return 0;
    }
    sact=nact+opt; eact=nact-1;
  }

  /* Enter loop over all actions */
  for (i=sact; i<=eact; i++) {

    /***** HERE - Adrian's code below to skip orders which, for some
	   reason, no longer have enough valid data *****/
    /***** Update in Version 0.60, 20/07/2012: This code is the root
	   cause of some segmentation faults I've been getting when
	   running existing UPL files. Removing it for now. ******/
    /* Detect zeroed-out orders before performing any further actions */
    /*
    if (spec[act[i].i[0]].or[act[i].i[1]].nuse==0) {
      nferrormsg("UVES_past_actions(): Action %d: No valid pixels in order %d in file %d,\n\t%s\n",i+1,
		 act[i].i[1]+1,act[i].i[0]+1,spec[act[i].i[0]].file);
      act[i].val = 0; // ??? can I do this?
      continue; // and this?
    }
    */

    if (act[i].val) {
      /* Some actions require finding pixels on the combined spectrum ... */
      if (act[i].act<SOACT) {
	if ((cp1=idxdval(cspec->rwl,cspec->np,act[i].d[0]))==-1)
	  errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[0]);
	if ((cp2=idxdval(&(cspec->rwl[cp1]),cspec->np-cp1,act[i].d[1]))==-1) {
	  if (act[i].d[1]>cspec->rwl[cspec->np-1]) cp2=cspec->np-1-cp1;
	  else
	    errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[1]);
	}
	cp2+=cp1;
	/* Some actions require finding pixels on individual orders as well */
	if (act[i].act==COACT || act[i].act==UOACT || act[i].act==FOACT ||
	    act[i].act==IOACT) {
	  if (cp2>=spec[act[i].i[0]].or[act[i].i[1]].csidx &&
	      cp1<=spec[act[i].i[0]].or[act[i].i[1]].ceidx) {
	    cp1=MAX(cp1,spec[act[i].i[0]].or[act[i].i[1]].csidx);
	    cp2=MIN(cp2,spec[act[i].i[0]].or[act[i].i[1]].ceidx);
	  } else {
	    nferrormsg("UVES_past_actions(): Action %d: Wavelength range\n\
\tof clip-window, %lf-%lf, does not overlap with that\n\
\tof order %d in file %d,\n\t%s,\n\twith wavelength range %lf - %lf. Skipping this action.\n",
		       i+1,act[i].d[0],act[i].d[1],act[i].i[1]+1,act[i].i[0]+1,spec[act[i].i[0]].file,
		       cspec->wl[spec[act[i].i[0]].or[act[i].i[1]].csidx],
		       cspec->wl[spec[act[i].i[0]].or[act[i].i[1]].ceidx]);
	    continue;
	  }
	  op1=cp1-spec[act[i].i[0]].or[act[i].i[1]].csidx;
	  op2=cp2-spec[act[i].i[0]].or[act[i].i[1]].csidx;
	}
      }
      
      /* Enter decision tree for different actions */
      switch (act[i].act) {
      case COACT:
	/** Clip pixel from order **/
	/* Flag those pixels to be cosmetically clipped */
	for (j=op1; j<=op2; j++) {
	  if (spec[act[i].i[0]].or[act[i].i[1]].rdfl[j]>act[i].d[2] &&
	      spec[act[i].i[0]].or[act[i].i[1]].rdfl[j]<act[i].d[3]) {
	    spec[act[i].i[0]].or[act[i].i[1]].ordst[j]=
	      spec[act[i].i[0]].or[act[i].i[1]].rdst[j];
	    spec[act[i].i[0]].or[act[i].i[1]].rdst[j]=CCLIP;
	  }
	}
	break;
      case CCACT:
	/** Clip pixel from combined spectrum **/
	/* Flag those pixels to be cosmetically clipped */
	for (j=cp1; j<=cp2; j++) {
	  if (cspec->fl[j]>act[i].d[2] && cspec->fl[j]<act[i].d[3]) {
	    cspec->ost[j]=cspec->st[j]; cspec->st[j]=CCLIP;
	    cspec->no[j]=1.0; cspec->ne[j]=cspec->nf[j]=-INFIN;
	  }
	}
	break;
      case UOACT:
	/** Unclip pixel from order **/
	/* Flag those pixels to be cosmetically unclipped */
	for (j=op1; j<=op2; j++) {
	  if (spec[act[i].i[0]].or[act[i].i[1]].rdfl[j]>act[i].d[2] &&
	      spec[act[i].i[0]].or[act[i].i[1]].rdfl[j]<act[i].d[3]) {
	    spec[act[i].i[0]].or[act[i].i[1]].ordst[j]=
	      spec[act[i].i[0]].or[act[i].i[1]].rdst[j];
	    spec[act[i].i[0]].or[act[i].i[1]].rdst[j]=1;
	  }
	}
	break;
      case UCACT:
	/** Unclip pixel from combined spectrum **/
	/* Flag those pixels to be cosmetically unclipped */
	for (j=cp1; j<=cp2; j++) {
	  if (cspec->fl[j]>act[i].d[2] && cspec->fl[j]<act[i].d[3]) {
	    cspec->ost[j]=cspec->st[j]; cspec->st[j]=1;
	    cspec->no[j]=cspec->fl[j]/cspec->co[j];
	    cspec->ne[j]=cspec->er[j]/cspec->co[j];
	    cspec->nf[j]=cspec->ef[j]/cspec->co[j];
	  }
	}
	break;
      case FOACT:
	/** Fit new continuum to single order **/
	/* Find correct pixels, first in combined spectrum, then in
	   specific order */
	if ((cp3=idxdval(&(cspec->rwl[cp1]),cspec->np-cp1,act[i].d[2]))==-1) {
	  if (act[i].d[2]>cspec->rwl[cspec->np-1]) cp3=cspec->np-2-cp1;
	  else errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[2]);
	}
	cp3+=cp1;
	if ((cp4=idxdval(&(cspec->rwl[cp3]),cspec->np-cp3,act[i].d[3]))==-1) {
	  if (act[i].d[3]>cspec->rwl[cspec->np-2]) cp4=cspec->np-2-cp3;
	  errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[3]);
	}
	cp4+=cp3;
	cp3=MAX(cp3,(spec[act[i].i[0]].or[act[i].i[1]].csidx+1));
	cp4=MIN(cp4,(spec[act[i].i[0]].or[act[i].i[1]].ceidx-1));
	op3=cp3-spec[act[i].i[0]].or[act[i].i[1]].csidx;
	op4=cp4-spec[act[i].i[0]].or[act[i].i[1]].csidx;
	op2=(MIN(op2,(spec[act[i].i[0]].or[act[i].i[1]].nrdp-1)));
	op4=(MIN(op4,(spec[act[i].i[0]].or[act[i].i[1]].nrdp-2)));
	/* Allocate memory for the temporary fitting arrays */
	fnp=op2-op1+1;
	if ((fco=darray(fnp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fco array of size %d",fnp);
	if ((fst=iarray(fnp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fst array of size %d",fnp);
	/* Fill status array with initial values */
	for (j=0,k=op1; j<fnp; j++,k++)
	  fst[j]=(spec[act[i].i[0]].or[act[i].i[1]].rdst[k]==SCLIP ||
		  spec[act[i].i[0]].or[act[i].i[1]].rdst[k]==LCLIP) ? 1 :
	    spec[act[i].i[0]].or[act[i].i[1]].rdst[k];	
	/* Now include information from deselected or reselected regions */
	for (j=0; j<act[i].nxyp; j++) {
	  /* Find correct pixels, first in combined spectrum, then in order */
	  if ((fp1=idxdval(cspec->rwl,cspec->np,act[i].xyp[j].x1))==-1)
	    errormsg("UVES_past_actions(): Action %d: XYP %d: Cannot find\n\
\tpixel near wavelength %lf in combined spectrum",i+1,j+1,act[i].xyp[j].x1);
	  if ((fp2=
	       idxdval(&(cspec->rwl[fp1]),cspec->np-fp1,act[i].xyp[j].x2))==-1) {
	    if (act[i].xyp[j].x2>cspec->rwl[cspec->np-1]) fp2=cspec->np-1;
	    else errormsg("UVES_past_actions(): Action %d: XYP %d: Cannot find\n\
\tpixel near wavelength %lf in combined spectrum",i+1,j+1,act[i].xyp[j].x2);
	  } else fp2+=fp1;
	  fp1-=spec[act[i].i[0]].or[act[i].i[1]].csidx+op1;
	  fp2-=spec[act[i].i[0]].or[act[i].i[1]].csidx+op1;
	  /* Fill with information from (de/re)selected regions */
	  if (fp2>=0 && fp1<fnp) {
	    fp1=MAX(fp1,0); fp2=MIN(fp2,fnp-1);
	    for (k=fp1,l=fp1+op1; k<=fp2; k++,l++)
	      if (spec[act[i].i[0]].or[act[i].i[1]].rdfl[l]>act[i].xyp[j].y1 &&
		  spec[act[i].i[0]].or[act[i].i[1]].rdfl[l]<act[i].xyp[j].y2)
		fst[k]=(act[i].xyp[j].i) ? 1 : PCLIP;
	  }
	}
	/* Do a fit for the selected data */
	frs=UVES_confit(&(spec[act[i].i[0]].or[act[i].i[1]].rdfl[op1]),
			&(spec[act[i].i[0]].or[act[i].i[1]].rder[op1]),fst,fnp,
			act[i].i[2],act[i].i[3],act[i].d[4],act[i].d[5],pctl,0,fco);
	if (frs==-1 || frs==0) {
	  warnmsg("UVES_past_actions(): Action %d: New continuum fit failed",i+1);
	  break;
	}
	/* Blend two continua together at ends over the blend region */
	fbsnp=op3-op1+1; fbenp=op2-op4+1;
	for (j=0,k=op1; j<fbsnp; j++,k++) {
	  frac=(double)(j+1)/(double)fbsnp;
	  fco[j]=frac*fco[j]+(1.0-frac)*spec[act[i].i[0]].or[act[i].i[1]].rdco[k];
	}
	for (j=0,k=fnp-1,l=op2; j<fbenp; j++,k--,l--) {
	  frac=(double)(j+1)/(double)fbenp;
	  fco[k]=frac*fco[k]+(1.0-frac)*spec[act[i].i[0]].or[act[i].i[1]].rdco[l];
	}
	/* Record the order pixels to reflect new continuum fit */
	if (par->thar<=1) {
	  for (j=0,k=op1; j<fnp; j++,k++) {
	    scale=spec[act[i].i[0]].or[act[i].i[1]].rdco[k]/fco[j];
	    spec[act[i].i[0]].or[act[i].i[1]].ordco[k]=scale;
	    spec[act[i].i[0]].or[act[i].i[1]].rdco[k]=fco[j];
	    spec[act[i].i[0]].or[act[i].i[1]].rdfl[k]*=scale;
	    spec[act[i].i[0]].or[act[i].i[1]].rder[k]*=scale;
	    spec[act[i].i[0]].or[act[i].i[1]].rdef[k]*=scale;
	    spec[act[i].i[0]].or[act[i].i[1]].rdme[k]*=scale;
	    /* Scale alternative arrays for pre-v0.74 backwards compatibility */
	    if (par->version<0.74 && spec[act[i].i[0]].or[act[i].i[1]].rdfl_a!=NULL) {
	      spec[act[i].i[0]].or[act[i].i[1]].rdfl_a[k]*=scale;
	      spec[act[i].i[0]].or[act[i].i[1]].rder_a[k]*=scale;
	      spec[act[i].i[0]].or[act[i].i[1]].rdef_a[k]*=scale;
	      spec[act[i].i[0]].or[act[i].i[1]].rdme_a[k]*=scale;
	    }
	  }
	}
	else {
	  for (j=0,k=op1; j<fnp; j++,k++) {
	    spec[act[i].i[0]].or[act[i].i[1]].ordco[k]=
	      spec[act[i].i[0]].or[act[i].i[1]].rdco[k];
	    spec[act[i].i[0]].or[act[i].i[1]].rdco[k]=fco[j];
	  }
	}
	/* Free temporary fitting arrays*/
	free(fco); free(fst);
	break;
      case FCACT:
	/** Fit new continuum to section of combined spectrum **/
	/* Find correct pixels */
	if ((cp3=idxdval(&(cspec->rwl[cp1]),cspec->np-cp1,act[i].d[2]))==-1)
	  errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[2]);
	cp3+=cp1;
	if ((cp4=idxdval(&(cspec->rwl[cp3]),cspec->np-cp3,act[i].d[3]))==-1) {
	  if (act[i].d[3]>cspec->rwl[cspec->np-2]) cp4=cspec->np-2-cp3;
	  else errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[3]);
	}
	cp4+=cp3;
	/* Allocate memory for the temporary fitting arrays */
	fnp=cp2-cp1+1;
	if ((fco=darray(fnp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fco array of size %d",fnp);
	if ((fst=iarray(fnp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fst array of size %d",fnp);
	/* Fill status array with initial values */
	for (j=0,k=cp1; j<fnp; j++,k++) fst[j]=cspec->st[k];
	/* Now include information from deselected or reselected regions */
	for (j=0; j<act[i].nxyp; j++) {
	  /* Find correct pixels, first in combined spectrum, then in order */
	  if ((fp1=idxdval(cspec->rwl,cspec->np,act[i].xyp[j].x1))==-1)
	    errormsg("UVES_past_actions(): Action %d: XYP %d: Cannot find\n\
\tpixel near wavelength %lf in combined spectrum",i+1,j+1,act[i].xyp[j].x1);
	  if ((fp2=idxdval(&(cspec->rwl[fp1]),cspec->np-fp1,act[i].xyp[j].x2))
	      ==-1) {
	    if (act[i].xyp[j].x2>cspec->rwl[cspec->np-1]) fp2=cspec->np-1;
	    else errormsg("UVES_past_actions(): Action %d: XYP %d: Cannot find\n\
\tpixel near wavelength %lf in combined spectrum",i+1,j+1,act[i].xyp[j].x2);
	  } else fp2+=fp1;
	  fp1-=cp1; fp2-=cp1;
	  /* Fill with information from (de/re)selected regions */
	  if (fp2>=0 && fp1<fnp) {
	    fp1=MAX(fp1,0); fp2=MIN(fp2,fnp-1);
	    for (k=fp1,l=fp1+cp1; k<=fp2; k++,l++)
	      if (cspec->fl[l]>act[i].xyp[j].y1 && cspec->fl[l]<act[i].xyp[j].y2)
		fst[k]=(act[i].xyp[j].i) ? 1 : PCLIP;
	  }
	}
	/* Do a fit for the selected data */
	frs=UVES_confit(&(cspec->fl[cp1]),&(cspec->er[cp1]),fst,fnp,act[i].i[0],
			act[i].i[1],act[i].d[4],act[i].d[5],pctl,0,fco);
	if (frs==-1 || frs==0) {
	  warnmsg("UVES_past_actions(): Action %d: New continuum fit failed",i+1);
	  break;
	}
	/* Blend two continua together at ends over the blend region */
	fbsnp=cp3-cp1+1; fbenp=cp2-cp4+1;
	for (j=0,k=cp1; j<fbsnp; j++,k++) {
	  frac=(double)(j+1)/(double)fbsnp;
	  fco[j]=frac*fco[j]+(1.0-frac)*cspec->co[k];
	}
	for (j=0,k=fnp-1,l=cp2; j<fbenp; j++,k--,l--) {
	  frac=(double)(j+1)/(double)fbenp;
	  fco[k]=frac*fco[k]+(1.0-frac)*cspec->co[l];
	}
	/* Find which orders contribute to the combined spectrum over the
	   relevant wavelength range */
	for (j=0,l=0; j<nspec; j++) {
	  for (k=0; k<spec[j].nor; k++) {
	    if (spec[j].or[k].nuse>=MINUSE &&
		cp2>=spec[j].or[k].csidx && cp1<=spec[j].or[k].ceidx) {
	      con[l][0]=j; con[l][1]=k; con[l][2]=MAX(cp1,spec[j].or[k].csidx);
	      con[l++][3]=MIN(cp2,spec[j].or[k].ceidx);
	    }
	  }
	}
	ncon=l;
	/* Record the combined spectrum and contributing order pixels to
	   reflect new continuum fit */
	if (par->thar<=1) {
	  for (j=0,k=cp1; j<fnp; j++,k++) {
	    cspec->oco[k]=cspec->co[k]; scale=cspec->co[k]/fco[j];
	    cspec->co[k]=fco[j]; cspec->no[k]*=scale;
	    cspec->ne[k]*=scale; cspec->nf[k]*=scale;
	    for (l=0; l<ncon; l++) {
	      if (con[l][2]<=k && con[l][3]>=k) {
		idx=k-spec[con[l][0]].or[con[l][1]].csidx;
		spec[con[l][0]].or[con[l][1]].rdco[idx]=cspec->co[k];
	      }
	    }
	  }
	}
	else {
	  for (j=0,k=cp1; j<fnp; j++,k++) {
	    cspec->oco[k]=cspec->co[k]; scale=cspec->co[k]-fco[j];
	    cspec->co[k]=fco[j]; cspec->no[k]+=scale;
	    for (l=0; l<ncon; l++) {
	      if (con[l][2]<=k && con[l][3]>=k) {
		idx=k-spec[con[l][0]].or[con[l][1]].csidx;
		spec[con[l][0]].or[con[l][1]].rdco[idx]=cspec->co[k];
	      }
	    }
	  }
	}
	/* Free temporary fitting arrays*/
	free(fco); free(fst);
	break;
      case IOACT:
	/** Interpolate new continuum to single order **/
	/* Find correct pixels, first in combined spectrum, then in
	   specific order */
	if ((cp3=idxdval(&(cspec->rwl[cp1]),cspec->np-cp1,act[i].d[2]))==-1) {
	  if (act[i].d[2]>cspec->rwl[cspec->np-1]) cp3=cspec->np-2-cp1;
	  else errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[2]);
	}
	cp3+=cp1;
	if ((cp4=idxdval(&(cspec->rwl[cp3]),cspec->np-cp3,act[i].d[3]))==-1) {
	  if (act[i].d[3]>cspec->rwl[cspec->np-2]) cp4=cspec->np-2-cp3;
	  errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[3]);
	}
	cp4+=cp3;
	cp3=MAX(cp3,(spec[act[i].i[0]].or[act[i].i[1]].csidx+1));
	cp4=MIN(cp4,(spec[act[i].i[0]].or[act[i].i[1]].ceidx-1));
	op3=cp3-spec[act[i].i[0]].or[act[i].i[1]].csidx;
	op4=cp4-spec[act[i].i[0]].or[act[i].i[1]].csidx;
	op2=(MIN(op2,(spec[act[i].i[0]].or[act[i].i[1]].nrdp-1)));
	op4=(MIN(op4,(spec[act[i].i[0]].or[act[i].i[1]].nrdp-2)));
	/* Allocate memory for the temporary fitting arrays */
	fnp=op2-op1+1;
	if ((fco=darray(fnp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fco array of size %d",fnp);
	if ((fx=darray(act[i].nxyp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fx array of size %d",act[i].nxyp);
	if ((fy=darray(act[i].nxyp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fy array of size %d",act[i].nxyp);
	if ((fy2=darray(act[i].nxyp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fy2 array of size %d",act[i].nxyp);
	/* Fill spline node arrays */
	for (j=0; j<act[i].nxyp; j++) {
	  fx[j]=act[i].xyp[j].x1; fy[j]=act[i].xyp[j].y1;
	}
	/* Do spline fit */
	if (!spline(fx,fy,act[i].nxyp,0.0,0.0,fy2,1))
	  errormsg("UVES_past_actions(): Action %d: Error calculating spline",i+1);
	for (j=0,k=cp1; j<fnp; j++,k++)
	  fco[j]=splint(fx,fy,fy2,act[i].nxyp,cspec->wl[k]);
	/* Blend two continua together at ends over the blend region */
	fbsnp=op3-op1+1; fbenp=op2-op4+1;
	for (j=0,k=op1; j<fbsnp; j++,k++) {
	  frac=(double)(j+1)/(double)fbsnp;
	  fco[j]=frac*fco[j]+(1.0-frac)*spec[act[i].i[0]].or[act[i].i[1]].rdco[k];
	}
	for (j=0,k=fnp-1,l=op2; j<fbenp; j++,k--,l--) {
	  frac=(double)(j+1)/(double)fbenp;
	  fco[k]=frac*fco[k]+(1.0-frac)*spec[act[i].i[0]].or[act[i].i[1]].rdco[l];
	}
	/* Record the order pixels to reflect new continuum fit */
	if (par->thar<=1) {
	  for (j=0,k=op1; j<fnp; j++,k++) {
	    scale=spec[act[i].i[0]].or[act[i].i[1]].rdco[k]/fco[j];
	    spec[act[i].i[0]].or[act[i].i[1]].ordco[k]=scale;
	    spec[act[i].i[0]].or[act[i].i[1]].rdco[k]=fco[j];
	    spec[act[i].i[0]].or[act[i].i[1]].rdfl[k]*=scale;
	    spec[act[i].i[0]].or[act[i].i[1]].rder[k]*=scale;
	    spec[act[i].i[0]].or[act[i].i[1]].rdef[k]*=scale;
	    spec[act[i].i[0]].or[act[i].i[1]].rdme[k]*=scale;
	    /* Scale alternative arrays for pre-v0.74 backwards compatibility */
	    if (par->version<0.74 && spec[act[i].i[0]].or[act[i].i[1]].rdfl_a!=NULL) {
	      spec[act[i].i[0]].or[act[i].i[1]].rdfl_a[k]*=scale;
	      spec[act[i].i[0]].or[act[i].i[1]].rder_a[k]*=scale;
	      spec[act[i].i[0]].or[act[i].i[1]].rdef_a[k]*=scale;
	      spec[act[i].i[0]].or[act[i].i[1]].rdme_a[k]*=scale;
	    }
	  }
	}
	else {
	  for (j=0,k=op1; j<fnp; j++,k++) {
	    spec[act[i].i[0]].or[act[i].i[1]].ordco[k]=
	      spec[act[i].i[0]].or[act[i].i[1]].rdco[k];
	    spec[act[i].i[0]].or[act[i].i[1]].rdco[k]=fco[j];
	  }
	}
	/* Free temporary fitting arrays*/
	free(fco); free(fx); free(fy); free(fy2);
	break;
      case ICACT:
	/** Fit new continuum to section of combined spectrum **/
	/* Find correct pixels */
	if ((cp3=idxdval(&(cspec->rwl[cp1]),cspec->np-cp1,act[i].d[2]))==-1)
	  errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[2]);
	cp3+=cp1;
	if ((cp4=idxdval(&(cspec->rwl[cp3]),cspec->np-cp3,act[i].d[3]))==-1) {
	  if (act[i].d[3]>cspec->rwl[cspec->np-2]) cp4=cspec->np-2-cp3;
	  else errormsg("UVES_past_actions(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",i+1,act[i].d[3]);
	}
	cp4+=cp3;
	/* Allocate memory for the temporary fitting arrays */
	fnp=cp2-cp1+1;
	if ((fco=darray(fnp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fco array of size %d",fnp);
	if ((fx=darray(act[i].nxyp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fx array of size %d",act[i].nxyp);
	if ((fy=darray(act[i].nxyp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fy array of size %d",act[i].nxyp);
	if ((fy2=darray(act[i].nxyp))==NULL)
	  errormsg("UVES_past_actions(): Could not allocate memory\n\
\tfor fy2 array of size %d",act[i].nxyp);
	/* Fill spline node arrays */
	for (j=0; j<act[i].nxyp; j++) {
	  fx[j]=act[i].xyp[j].x1; fy[j]=act[i].xyp[j].y1;
	}
	/* Do spline fit */
	if (!spline(fx,fy,act[i].nxyp,0.0,0.0,fy2,1))
	  errormsg("UVES_past_actions(): Action %d: Error calculating spline",i+1);
	for (j=0,k=cp1; j<fnp; j++,k++)
	  fco[j]=splint(fx,fy,fy2,act[i].nxyp,cspec->wl[k]);
	/* Blend two continua together at ends over the blend region */
	fbsnp=cp3-cp1+1; fbenp=cp2-cp4+1;
	for (j=0,k=cp1; j<fbsnp; j++,k++) {
	  frac=(double)(j+1)/(double)fbsnp;
	  fco[j]=frac*fco[j]+(1.0-frac)*cspec->co[k];
	}
	for (j=0,k=fnp-1,l=cp2; j<fbenp; j++,k--,l--) {
	  frac=(double)(j+1)/(double)fbenp;
	  fco[k]=frac*fco[k]+(1.0-frac)*cspec->co[l];
	}
	/* Find which orders contribute to the combined spectrum over the
	   relevant wavelength range */
	for (j=0,l=0; j<nspec; j++) {
	  for (k=0; k<spec[j].nor; k++) {
	    if (spec[j].or[k].nuse>=MINUSE &&
		cp2>=spec[j].or[k].csidx && cp1<=spec[j].or[k].ceidx) {
	      con[l][0]=j; con[l][1]=k;
	      con[l][2]=MAX(cp1,spec[j].or[k].csidx);
	      con[l++][3]=MIN(cp2,spec[j].or[k].ceidx);
	    }
	  }
	}
	ncon=l;
	/* Record the combined spectrum and contributing order pixels to
	   reflect new continuum fit */
	if (par->thar<=1) {
	  for (j=0,k=cp1; j<fnp; j++,k++) {
	    cspec->oco[k]=cspec->co[k]; scale=cspec->co[k]/fco[j];
	    cspec->co[k]=fco[j]; cspec->no[k]*=scale;
	    cspec->ne[k]*=scale; cspec->nf[k]*=scale;
	    for (l=0; l<ncon; l++) {
	      if (con[l][2]<=k && con[l][3]>=k) {
		idx=k-spec[con[l][0]].or[con[l][1]].csidx;
		spec[con[l][0]].or[con[l][1]].rdco[idx]=cspec->co[k];
	      }
	    }
	  }
	}
	else {
	  for (j=0,k=cp1; j<fnp; j++,k++) {
	    cspec->oco[k]=cspec->co[k]; scale=cspec->co[k]-fco[j];
	    cspec->co[k]=fco[j]; cspec->no[k]+=scale;
	    for (l=0; l<ncon; l++) {
	      if (con[l][2]<=k && con[l][3]>=k) {
		idx=k-spec[con[l][0]].or[con[l][1]].csidx;
		spec[con[l][0]].or[con[l][1]].rdco[idx]=cspec->co[k];
	      }
	    }
	  }
	}
	/* Free temporary fitting arrays*/
	free(fco); free(fx); free(fy); free(fy2);
	break;
      case SOACT:
	/** Scale order flux, error and continuum arrays **/
	for (j=0; j<spec[act[i].i[0]].or[act[i].i[1]].nrdp; j++) {
	  spec[act[i].i[0]].or[act[i].i[1]].rdfl[j]*=act[i].d[0];
	  spec[act[i].i[0]].or[act[i].i[1]].rder[j]*=act[i].d[0];
	  spec[act[i].i[0]].or[act[i].i[1]].rdef[j]*=act[i].d[0];
	  spec[act[i].i[0]].or[act[i].i[1]].rdme[j]*=act[i].d[0];
	}
	/* Scale alternative arrays for pre-v0.74 backwards compatibility */
	if (par->version<0.74 && spec[act[i].i[0]].or[act[i].i[1]].rdfl_a!=NULL) {
	  for (j=0; j<spec[act[i].i[0]].or[act[i].i[1]].nrdp; j++) {
	    spec[act[i].i[0]].or[act[i].i[1]].rdfl_a[j]*=act[i].d[0];
	    spec[act[i].i[0]].or[act[i].i[1]].rder_a[j]*=act[i].d[0];
	    spec[act[i].i[0]].or[act[i].i[1]].rdef_a[j]*=act[i].d[0];
	    spec[act[i].i[0]].or[act[i].i[1]].rdme_a[j]*=act[i].d[0];
	  }
	}
	spec[act[i].i[0]].or[act[i].i[1]].oscl=spec[act[i].i[0]].or[act[i].i[1]].scl;
	spec[act[i].i[0]].or[act[i].i[1]].scl*=act[i].d[0];
	break;
      case ARACT:
	/** Auto-rescale orders after manual scaling **/
	status=UVES_rescale_region(spec,nspec,cspec,act,i,act[i].d[0],act[i].d[1],
				   par);
	switch (status) {
	case 0:
	  errormsg("UVES_past_actions(): Error returned from\n\
\tUVES_rescale_region()");
	  break;
	case 2: warnmsg("UVES_past_actions(): Action %d: No manually\n\
\trescaled orders after last recombination",i+1); break;
	}
	break;
      case NCACT:
	/** Create a new continuumm or entire spectrum based on
	    command-line parameters, including combmeth **/
	if (!par->combmeth) {
	  if (!UVES_cspec_cont(spec,nspec,cspec,par)) {
	    nferrormsg("UVES_past_actions(): Error returned from\n\
\tUVES_cspec_cont() when attempting to reset continuum"); return 0;
	  }
	} else {
	  /* Reset continuum before rederiving it */
	  for (j=0; j<nspec; j++) {
	    for (k=0; k<spec[j].nor; k++) {
	      if (spec[j].or[k].nuse>=MINUSE) {
		for (l=0; l<spec[j].or[k].nrdp; l++)
		  if (spec[j].or[k].rdst[l]==LCLIP && spec[j].or[k].rdst[l]==SCLIP)
		    spec[j].or[k].rdst[l]=1;
	      }
	    }
	  }
	  /* Make a rough continuum for each order of every spectrum */
	  if (!UVES_order_cont(spec,nspec,par))
	    errormsg("UVES_past_actions(): Unknown error returned from\n\
\tUVES_order_cont()");
	  /* Make a unified continuum for the combined spectrum */
	  if (!UVES_combine_cont(spec,nspec,cspec,1,par))
	    errormsg("UVES_past_actions(): Unknown error returned from\n\
\tUVES_combine_cont()");
	}
	break;
      case AAACT:
	/** Switch alternative arrays to main arrays. Do this by
	    freeing memory for main arrays, changing their pointers to
	    the alternative arrays. Then set the alternative flux
	    array pointers to NULL to signal to other subroutines to
	    not keep them up to date anymore **/
	if (par->version<0.74) {
	  for (j=0; j<nspec; j++) {
	    for (k=0; k<spec[j].nor; k++) {
	      if (spec[j].or[k].nuse>=MINUSE && spec[j].or[k].rdfl_a!=NULL) {
		free(spec[j].or[k].rdfl); spec[j].or[k].rdfl=spec[j].or[k].rdfl_a;
		spec[j].or[k].rdfl_a=NULL;
		free(spec[j].or[k].rder); spec[j].or[k].rder=spec[j].or[k].rder_a;
		spec[j].or[k].rder_a=NULL;
		free(spec[j].or[k].rdef); spec[j].or[k].rdef=spec[j].or[k].rdef_a;
		spec[j].or[k].rdef_a=NULL;
		free(spec[j].or[k].rdme); spec[j].or[k].rdme=spec[j].or[k].rdme_a;
		spec[j].or[k].rdme_a=NULL;
	} } } }
	break;
	/* Exit decision tree for different action types */
      }

      /* Recombine the spectra if required */
      if (act[i].rcmb && act[i].act!=ARACT)
	if (!UVES_combine_spec(spec,nspec,cspec,par))
	  errormsg("UVES_past_actions(): Action %d: Problem combining\n\
\tspectra in UVES_combine_spec()",i+1);

      /* Register a single action in one go for each action */
      act[i].nordact=1;
    }
    /* Exit loop over all actions */
  }

  /* Clean up */
  free(*con); free(con);

  return 1;

}
