/****************************************************************************
* Plot out final combined spectrum
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "UVES_popler.h"
#include "fit.h"
#include "stats.h"
#include "sort.h"
#include "gamm.h"
#include "input.h"
#include "memory.h"
#include "const.h"
#include "error.h"

int UVES_plot_cspec(spectrum *spec, int nspec, cspectrum *cspec, cplot *cp,
		    rplot *rp, action **act, int *nact, int *nact_save,
		    params *par) {

  double         frac=0.0,rejsigu=0.0,rejsigl=0.0,scale=0.0,sum=0.0,pe=0.0;
  double         ew=0.0,eew=0.0,cwl=0.0,cwin=0.0;
  double         wl1=0.0,wl2=0.0;
  double         sscale=SSCALE,bscale=BSCALE;
  double         *dat=NULL,*err=NULL,*efl=NULL,*med=NULL;
  double         *fx=NULL,*fy=NULL,*fy2=NULL;
  double         *fco=NULL;
  float          swl=0.0,ewl=0.0;
  float          flmin=0.0,flmax=0.0;
  float          csqmin=0.0,csqmax=0.0;
  float          ncbmin=0.0,ncbmax=0.0;
  float          shflmin=INFIN,shflmax=0.0;
  float          pgx=0.0,pgy=0.0,pgxo=0.0,pgyo=0.0;
  float          step=0.0,clickx=0.0,fscale=0.0,cdist=0.0,ftemp=0.0;
  float          ftmp;
  float          *shwl=NULL,*shfl=NULL;
  float          *wl=NULL,*fl=NULL,*er=NULL,*co=NULL,*csq=NULL,*ccsq=NULL;
  float          *no=NULL,*ne=NULL,*ncb=NULL,*nccb=NULL;
  float          *or=NULL;
  int            nxyp=0,cact=0,ndat=0;
  int            sp=0,ep=0,cnp=0,np=0,shnp=0;
  int            csp=0,cep=0,cosp,coep=0;
  int            nor=0,ncon=0,npcon=0,nscon=0,tscon=0,idx=0;
  int            fsp=0,fep=0,fnp=0,fbsp=0,fbep=0,fbsnp=0,fbenp=0,fosp=0,foep=0;
  int            plval=0,status=0,soact=0,sclinit=-1,con_redo=1;
  int            fit_typ=0,fit_ord=0,fit_rs=0,clicknp=0;
  int            i=0,j=0,k=0,l=0,m=0,firstfit=0,selpix=0;
  int            *fst=NULL;
  int            **con=NULL,**scon=NULL;
  char           pgch[SHORTLEN]="p";
  char           query[QUERYLEN]="\0";
  statset        stdat;
  twoxyp         *xyp=NULL;
  dbletwoint     *rank=NULL;
  extern plotenv plenv;

  /* Determine some global initial values and allocate memory for plotting
     arrays */
  cp->scalclip=par->scalclip; cp->scalerr=par->scalerr;
  cp->exit=0;
  np=cspec->np; shnp=(double)np/cp->sfac;
  cact=(par->replay) ? rp->cact : *nact-1;

  /* Allocate memory and fill arrays for low-dispersion version of spectrum */
  if ((shwl=farray(shnp))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor shwl plotting array of size %d",shnp);
  if ((shfl=farray(shnp))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor shfl plotting array of size %d",shnp);
  for (i=0,j=0; i<shnp; i++,j+=(int)cp->sfac) {
    shwl[i]=cspec->wl[j]; shfl[i]=cspec->fl[j];
  }
  for (i=0; i<np; i++) {
    if (cspec->st[i]>0) {
      shflmin=(MIN(cspec->fl[i],shflmin)); shflmax=(MAX(cspec->fl[i],shflmax));
    }
  }

  /* Allocate memory for plotting arrays */
  if ((wl=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor wl plotting array of size %d",np);
  if ((fl=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor fl plotting array of size %d",np);
  if ((er=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor er plotting array of size %d",np);
  if ((co=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor co plotting array of size %d",np);
  if ((no=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor no plotting array of size %d",np);
  if ((ne=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor ne plotting array of size %d",np);
  if ((csq=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor csq plotting array of size %d",np);
  if ((ccsq=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor ccsq plotting array of size %d",np);
  if ((ncb=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor ncb plotting array of size %d",np);
  if ((nccb=farray(np))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor nccb plotting array of size %d",np);
  /* Fill plotting arrays */
  for (i=0; i<np; i++) {
    wl[i]=cspec->wl[i]; fl[i]=cspec->fl[i]; er[i]=cspec->er[i];
    co[i]=cspec->co[i]; csq[i]=cspec->csq[i]; ccsq[i]=cspec->ccsq[i];
    no[i]=cspec->no[i]; ne[i]=cspec->ne[i];
    ncb[i]=(float)cspec->ncb[i]; nccb[i]=(float)cspec->nccb[i];
  }

  /* Set up plotting parameters */
  if (!cp->ep) {
    /** Previous combined plot info not set so set it up here **/
    /* Set wavelength and plotting limits */
    sp=MAX(0,(MIN(((np-cp->np)/2),(np-cp->np)))); ep=MIN((np-1),(sp+cp->np-1));
    cnp=ep-sp+1; swl=wl[sp]; ewl=wl[ep];
    flmax=1.1*fjmax(&(fl[sp]),cnp); flmin=MIN((-0.1*flmax),(fjmin(&(fl[sp]),cnp)));
  }
  else {
    sp=cp->sp; ep=cp->ep; swl=cp->swl; ewl=cp->ewl; cnp=cp->np;
    flmin=cp->ymn; flmax=cp->ymx;
  }
  plval=cp->plval;
  csqmin=fjmin(&(csq[sp]),cnp); csqmax=fjmax(&(csq[sp]),cnp);
  ncbmin=fjmin(&(ncb[sp]),cnp); ncbmax=fjmax(&(ncb[sp]),cnp);
  plenv.xmin[0]=shwl[0]; plenv.xmax[0]=shwl[shnp-1];
  plenv.ymax[0]=1.1*shflmax; plenv.ymin[0]=-0.1*plenv.ymax[0];
  plenv.xmin[1]=swl; plenv.xmax[1]=ewl;
  plenv.ymax[1]=1.1*(MAX(1.0,ncbmax)); plenv.ymin[1]=-0.1*plenv.ymax[1];
  plenv.ymax[2]=1.1*(MAX(0.5,(MIN(5.0,csqmax))));
  plenv.ymin[2]=-0.1*plenv.ymax[2];
  plenv.ymin[3]=flmin; plenv.ymax[3]=flmax;
  plenv.ymin[4]=-0.2; plenv.ymax[4]=1.3;

  /* Find longest order and allocate memory to order array */
  for (i=0,k=0; i<nspec; i++) {
    for (j=0; j<spec[i].nor; j++) {
      if (spec[i].or[j].nuse>=MINUSE) {
	if (spec[i].or[j].nrdp>nor) nor=spec[i].or[j].nrdp;
	k++;
      }
    }
  }
  if ((or=farray(nor))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor or plotting array of size %d",nor);

  /* Allocate memory to contribution matrix */
  if ((con=imatrix(k,7))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor con matrix of size %dx%d",k,7);

  /* Allocate memory to selected contribution matrix */
  if ((scon=imatrix(k,2))==NULL)
     errormsg("UVES_plot_cspec(): Cannot allocate memory\n\
\tfor scon matrix of size %dx%d",k,2);

  /* Allocate memory for rank array of structures */
  if (!(rank=(dbletwoint *)malloc((size_t)(k*sizeof(dbletwoint)))))
    errormsg("UVES_plot_cspec(): Cannot allocate memory for rank\n\
\tarray of length %d",k);

  /* Select correct PGPLOT device and advance to next page */
  cpgslct(cp->pgid); cpgpage();

  /** Enter plotting loop **/
  while (!cp->exit) {

    /* Begin plotting buffer */
    cpgbbuf(); cpgeras();

    /* Only allow plotting of contributing orders when some exist */
    if (!nspec) plval=0;

    /* Detect any changes in size of plotting window */
    cpgqvsz(1,&ftmp,&(plenv.wwidth),&ftmp,&(plenv.wasp)); plenv.wasp/=plenv.wwidth;

    /* Find which orders contribute to the combined spectrum over the
       relevant wavelength range */
    if (con_redo) {
      for (i=0,k=0; i<nspec; i++) {
	for (j=0; j<spec[i].nor; j++) {
	  if ((ep>=spec[i].or[j].csidx && sp<=spec[i].or[j].ceidx) &&
	      spec[i].or[j].nuse>=MINUSE) {
	    for (l=(MAX(sp,spec[i].or[j].csidx))-spec[i].or[j].csidx,m=0,sum=0.0;
		 l<=(MIN(ep,spec[i].or[j].ceidx))-spec[i].or[j].csidx; l++) {
	      if (spec[i].or[j].rdst[l]==1) {
		sum+=spec[i].or[j].rdfl[l]/spec[i].or[j].rder[l]; m++;
	      }
	    }
	    if (!cp->irank) {
	      rank[k].a=(m) ? sum/(double)m+DRNDTOL*(double)k : DRNDTOL*(double)k;
	    } else if (cp->irank==1) rank[k].a=(double)spec[i].or[j].crank;
	    else if (cp->irank==2) rank[k].a=(double)spec[i].or[j].ceidx;
	    rank[k].i=i; rank[k].j=j; k++;
	  }
	}
      }
      ncon=k;

      /* Sort the rank array so that we can rank plotted orders in
	 terms of SNR or combination rank */
      if (nspec) qsort(rank,ncon,sizeof(dbletwoint),qsort_dbletwointarray);

      /* Now rank contributing orders in terms of SNR or combination rank */
      for (i=0,j=1,soact=0; i<ncon; i++) {
	con[i][0]=rank[i].i; con[i][1]=rank[i].j;
	con[i][2]=MAX(sp,spec[con[i][0]].or[con[i][1]].csidx);
	con[i][3]=MIN(ep,spec[con[i][0]].or[con[i][1]].ceidx);
	con[i][4]=con[i][3]-con[i][2]+1;
	con[i][5]=j=(PMOD(j,1,14,2));
	if (tscon) {
	  for (k=0; k<nscon; k++) {
	    if (scon[k][0]==con[i][0] && scon[k][1]==con[i][1]) {
	      con[i][6]=1; break;
	    }
	  }
	  if (k==nscon) con[i][6]=0;
	} else if (par->replay && rp->exit>=3) {
	  if (con[i][0]==(*act)[rp->cact].i[0] &&
	      con[i][1]==(*act)[rp->cact].i[1]) {
	    if ((*act)[rp->cact].act==COACT || (*act)[rp->cact].act==UOACT ||
		(*act)[rp->cact].act==FOACT || (*act)[rp->cact].act==IOACT)
	      con[i][6]=plval=1;
	    else if ((*act)[rp->cact].act==SOACT) con[i][6]=soact=plval=1;
	  } else if (soact) con[i][6]=1;
	} else con[i][6]=1;
      }
    }

    /* Plot low-dispersion version of spectrum */
    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd0,cp->vpu0); cpgsci(0);
    cpgswin(plenv.xmin[0],plenv.xmax[0],plenv.ymin[0],plenv.ymax[0]);
    cpgsci(15);
    cpgrect(plenv.xmin[0],plenv.xmax[0],plenv.ymin[0],plenv.ymax[0]);
    cpgsci(5);
    cpgsch(0.3*plenv.ch); cpgbox("BCTS",0.0,0,"BC",0.0,0); cpgsch(plenv.ch);
    cpgsci(7); cpgsch(0.8*plenv.ch); cpgbox("N",0.0,0,"",0.0,0);
    cpgsci(3); cpgmtxt("B",2.5,0.5,0.5,plenv.xlab[0]);
    cpgsci(1); cpgsch(plenv.ch); cpgbin(shnp,shwl,shfl,1);
    cpgsci(3); cpgslw(2.0*plenv.lw); cpgsfs(2);
    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
    cpgslw(plenv.lw); cpgsfs(1);

    /* Plot the pixel-number spectra (before and after sigma-clipping) */
    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd1,cp->vpu1); cpgsci(0);
    cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[1],plenv.ymax[1]);
    cpgsci(15);
    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[1],plenv.ymax[1]);
    cpgsci(5); cpgsch(0.7*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
    cpgsci(7); cpgsch(0.8*plenv.ch); cpgbox("N",0.0,0,"",0.0,0);
    cpgsch(0.6*plenv.ch); cpgbox("",0.0,0,"NV",0.0,0);
    cpgsci(3); cpgsch(0.8*plenv.ch); cpglab(" ",plenv.ylab[1]," ");
    cpgsci(2); cpgbin(cnp,&(wl[sp]),&(ncb[sp]),1);
    cpgsci(1); cpgbin(cnp,&(wl[sp]),&(nccb[sp]),1);
    
    /* Plot the chisq. spectra (before and after sigma-clipping) */
    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd2,cp->vpu2); cpgsci(0);
    cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[2],plenv.ymax[2]);
    cpgsci(15);
    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[2],plenv.ymax[2]);
    cpgsci(5); cpgsch(0.7*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
    cpgsci(7); cpgsch(0.6*plenv.ch); cpgbox("",0.0,0,"NV",0.0,0);
    cpgsci(3); cpgsch(0.8*plenv.ch); cpglab(" ",plenv.ylab[2]," ");
    cpgsci(2); cpgbin(cnp,&(wl[sp]),&(csq[sp]),1);
    cpgsci(1); cpgbin(cnp,&(wl[sp]),&(ccsq[sp]),1);

    /* Plot the combined spectrum or contributing spectra */
    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd3,cp->vpu3); cpgsci(0);
    cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
    cpgsci(15);
    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
    cpgsci(5);
    cpgsch(0.9*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0); cpgsch(plenv.ch);
    cpgsci(3); cpgsch(0.8*plenv.ch); cpglab(" ",plenv.ylab[3]," ");
    if (plval) {
      npcon=0; for (i=0,l=0; i<ncon; i++) {
	if (con[i][6]) {
	  npcon++; l=i;
	  for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
	       j<con[i][4]; j++,k++)
	    or[j]=spec[con[i][0]].or[con[i][1]].rdfl[k];
	  cpgsci(con[i][5]); cpgbin(con[i][4],&(wl[con[i][2]]),or,1);
	  for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
	       j<con[i][4]; j++,k++) {
	    switch (spec[con[i][0]].or[con[i][1]].rdst[k]) {
	    case RCLIP: cpgpt(1,&(wl[con[i][2]+j]),&(or[j]),12); break;
	    case OCLIP: cpgpt(1,&(wl[con[i][2]+j]),&(or[j]),3); break;
	    case ACLIP: cpgpt(1,&(wl[con[i][2]+j]),&(or[j]),8); break;
	    case LCLIP: cpgpt(1,&(wl[con[i][2]+j]),&(or[j]),4); break;
	    case NCLIP: cpgpt(1,&(wl[con[i][2]+j]),&(or[j]),5); break;
	    case SCLIP: cpgpt(1,&(wl[con[i][2]+j]),&(or[j]),6); break;
	    case CCLIP: cpgpt(1,&(wl[con[i][2]+j]),&(or[j]),7); break;
	    }
	  }
	}
      }
      if (npcon==1) {
	for (j=0,k=con[l][2]-spec[con[l][0]].or[con[l][1]].csidx; j<con[l][4];
	     j++,k++) or[j]=spec[con[l][0]].or[con[l][1]].rdme[k];
	cpgsci((PMOD(con[l][5],3,14,2))); cpgline(con[l][4],&(wl[con[l][2]]),or);
	for (j=0,k=con[l][2]-spec[con[l][0]].or[con[l][1]].csidx; j<con[l][4];
	     j++,k++) or[j]=spec[con[l][0]].or[con[l][1]].rder[k];
	cpgsci((PMOD(con[l][5],1,14,2))); cpgline(con[l][4],&(wl[con[l][2]]),or);
	cpgsci((PMOD(con[l][5],2,14,2)));
      } else cpgsci(3);
      if (npcon) cpgline(cnp,&(wl[sp]),&(co[sp]));
    } else {
      cpgsci(1); cpgbin(cnp,&(wl[sp]),&(fl[sp]),1);
      cpgsci(2); cpgline(cnp,&(wl[sp]),&(er[sp]));
      cpgsci(3); cpgline(cnp,&(wl[sp]),&(co[sp])); cpgsci(7);
      for (i=sp; i<=ep; i++) {
	switch (cspec->st[i]) {
	case 1: break;
	case RCLIP: cpgpt(1,&(wl[i]),&(fl[i]),12); break;
	case NCLIP: cpgpt(1,&(wl[i]),&(fl[i]),5); break;
	case CCLIP: cpgpt(1,&(wl[i]),&(fl[i]),7); break;
	}
      }
    }
    cpgsci(5); cpgsls(4);
    cpgmove(plenv.xmin[1],0.0); cpgdraw(plenv.xmax[1],0.0); cpgsls(1);
    if (par->replay && rp->exit>=3 && (*act)[rp->cact].act<=UCACT) {
      cpgsci(7); cpgsfs(2);
      cpgrect((float)(*act)[rp->cact].d[0],(float)(*act)[rp->cact].d[1],
	      (float)(*act)[rp->cact].d[2],(float)(*act)[rp->cact].d[3]);
      cpgsfs(1);
    }
    
    /* Plot the normalized spectrum */
    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd4,cp->vpu4); cpgsci(0);
    cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[4],plenv.ymax[4]);
    cpgsci(15);
    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[4],plenv.ymax[4]);
    cpgsci(5); cpgsch(0.7*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
    cpgsci(7); cpgsch(0.6*plenv.ch); cpgbox("",0.0,0,"NV",0.0,0);
    cpgsci(3); cpgsch(0.8*plenv.ch); cpglab(" ",plenv.ylab[4]," ");
    cpgsci(1); cpgbin(cnp,&(wl[sp]),&(no[sp]),1);
    cpgsci(2); cpgline(cnp,&(wl[sp]),&(ne[sp]));
    cpgsci(5); cpgsch(plenv.ch); cpgsls(4);
    cpgmove(plenv.xmin[1],0.0); cpgdraw(plenv.xmax[1],0.0);
    cpgmove(plenv.xmin[1],1.0); cpgdraw(plenv.xmax[1],1.0); cpgsls(1);

    /* If required, plot a panel with buttons for switching on and off
       contributing spectra */
    if (plval && ncon>1) {
      cpgsvp(plenv.vpl,plenv.vpr,plenv.vpu,1.04*plenv.vpu);
      cpgswin(0.0,1.0,0.0,1.0); cpgslw(2.0*plenv.lw);
      cpgsch(plenv.ch); cpgsci(1);
      if (!cp->irank) cpgmtxt("B",-0.3,0.007,0.0,"S N R");
      else if (cp->irank==1) cpgmtxt("B",-0.3,0.013,0.0,"RANK");
      else if (cp->irank==2) cpgmtxt("B",-0.3,0.007,0.0,"O'LAP");
      step=0.9/(double)ncon;
      for (i=0; i<ncon; i++) {
	cpgsci(con[i][5]); if (!con[i][6]) cpgsfs(2);
	cpgrect(0.1+(float)i*step,0.1+(float)(i+1)*step,0.02,1.0);
	cpgsfs(1);
      }
      cpgslw(plenv.lw);
    }
    
    /* Flush the buffer */
    cpgebuf();

    /* No need to go further if replaying certain actions */
    if (par->replay && rp->exit==4) { rp->exit=0; break; }

    /** Main outer selection loop - waiting for navigation/editing
	instructions from user **/
    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd0,cp->vpu4);
    cpgswin(plenv.vpl,plenv.vpr,cp->vpd0,cp->vpu4);
    /* If in replay mode, decide whether interaction with Spectrum
       Navigator is necessary */
    if (par->replay && rp->exit==3) {
      if ((*act)[rp->cact].act==FOACT || (*act)[rp->cact].act==FCACT)
	sprintf(pgch,"f");
      else if ((*act)[rp->cact].act==IOACT || (*act)[rp->cact].act==ICACT)
	sprintf(pgch,"i");
    }
    /* If normal mode then await the user's instructions */
    else cpgband(0,0,0.0,0.0,&pgx,&pgy,pgch);
    /* Give usage instructions if asked */
    if (!strncmp(pgch,"?",1)) {
      if (plval) fprintf(stderr,"\
Options for upper selection panel:\n\
 Left mouse                : Selects single order to plot\n\
 Middle mouse              : Selects only higher S/N orders\n\
 Right mouse               : Toggles plotting of selected order\n\
 2 or 1                    : Down-scales order by factor %7.3lf or %7.2lf (act)\n\
 3 or 4                    : Up-scales order by factor %7.3lf or %7.2lf (act)\n\
 e                         : Selects orders in same exposure as selected order\n\
 i                         : Identify spectrum & order numbers & write to screen\n\
 l                         : Selects only lower S/N orders\n\
 n                         : Selects next highest S/N order & de-selects current;\n\
                                also works on main selection panel.\n\
 p                         : Selects next lowest S/N order & de-selects current;\n\
                                also works on main selection panel.\n\
 x                         : Selects orders in all spectra other than that selected\n",
			 sscale,bscale,sscale,bscale);
fprintf(stderr,"\
Options for main selection panel:\n\
 Left mouse on lower panel : Selects lower panel for navigation\n\
 Left mouse on flux panel  : Selects flux panel for navigation\n\
 Space bar                 : Refreshes the panel. Good after window resizing\n\
 N                         : Derive new automatic continuum for entire spectrum\n\
                               with command-line param.s depending on combmeth\n\
 a                         : Autorescale orders after manually scaling one or\n\
                               more (order scaling combination method only; act)\n\
 c                         : Adjust sig-clipping level for autorescaling orders\n\
 d                         : Delete or add spectra to plotted and saved comb.\n\
 e                         : Adjust scaling error threshold for autorescaling\n\
 f                         : Fit new continuum to combined spectrum or\n\
                               single order with normal polynomial fit (*; act)\n\
 i                         : Interpolate new continuum to combined spectrum or\n\
                               single order with spline function (*; act)\n\
 k                         : Keep selected spectra only when navigating\n\
 m                         : Mirror last action on currently displayed order(s);\n\
                               Works for order cont. fits and order scaling only.\n\
 o                         : Toggle ranking of orders in SNR or combination rank\n\
                               when using order scaling combination method\n\
 r                         : Recombine contributing spectra to form new\n\
                               combined spectrum (act)\n\
 s                         : Save actions in UPL file & save FITS 1-D spectrum\n\
 u                         : Undo last action\n\
 v                         : Plot source spectra with pixel status symbols\n\
 w                         : Toggle plotting of combined spectrum in white in\n\
                               front, behind or off when viewing source spectra\n\
 y                         : Re-scale y-axis to suit local (combined spectrum)\n\
                               extrema\n\
 .                         : Expand y-scale by factor of 2\n\
 ,                         : Shrink y-scale by factor of 2\n\
 >                         : Expand x-axis by factor of 2 (twice as many pixels)\n\
 <                         : Shrink x-axis by factor of 2 (half as many pixels)\n\
 q                         : Quit screen and exit program\n\n");
    }
    else if (!strncmp(pgch,"A",1) && pgy<cp->vpu0) {
      /* Use low-dispersion spectrum to navigate combined spectrum */
      cpgsvp(plenv.vpl,plenv.vpr,cp->vpd0,cp->vpu0);
      cpgswin(plenv.xmin[0],plenv.xmax[0],plenv.ymin[0],plenv.ymax[0]);
      cpgsci(7); cpgslw(2.0*plenv.lw);
      cpgsch(0.3*plenv.ch); cpgbox("BCTS",0.0,0,"BC",0.0,0); cpgsch(plenv.ch);
      cpgslw(plenv.lw); cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
      if (!strncmp(pgch,"?",1)) {
	fprintf(stderr,"\
Options for navigation using 'thumbnail' spectrum:\n\
 Left mouse  : Re-centre detailed spectrum (flux panel) at selected\n\
                 wavelength\n\
 Middle mouse: Re-centre detailed spectrum (as above) and re-size\n\
                 window to suit new local extrema\n\
 Right mouse : Select custom window region.\n\
                 Right mouse again selects second corner of new window,\n\
                 any other entry to abort selection.\n\n");
	cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
      }

      if (!strncmp(pgch,"D",1)) {
	pgx=(MIN(wl[np-1],(MAX(pgx,wl[0]))));
	sp=MAX(0,(MIN((idxfval(wl,np,pgx)-cp->np/2),(np-cp->np))));
	ep=sp+cp->np-1; cnp=ep-sp+1;
	plenv.xmin[1]=wl[sp]; plenv.xmax[1]=wl[ep];
	plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
	plenv.ymin[1]=-0.1*plenv.ymax[1];
	plenv.ymax[2]=1.1*(MAX(0.5,(MIN(5.0,(fjmax(&(csq[sp]),cnp))))));
	plenv.ymin[2]=-0.1*plenv.ymax[2];
	plenv.ymax[3]=fjmax(&(fl[sp]),cnp);
	plenv.ymin[3]=MIN((-0.1*plenv.ymax[3]),(fjmin(&(fl[sp]),cnp)));
	con_redo=1;
      }
      else if (!strncmp(pgch,"X",1)) {
	strcpy(pgch,"\0"); cpgsci(2); pgxo=pgx; pgyo=pgy;
	cpgband(2,0,pgxo,pgyo,&pgx,&pgy,pgch);
	if (!strncmp(pgch,"X",1)) {
	  if (pgxo>pgx) { FSWAP(pgxo,pgx); }
	  if (pgyo>pgy) { FSWAP(pgyo,pgy); }
	  pgxo=(MAX(pgxo,wl[0])); pgx=(MIN(pgx,wl[np-1]));
	  sp=MAX(0,(MIN((idxfval(wl,np,pgxo)),(np-2))));
	  ep=MAX(sp+1,(MIN((idxfval(wl,np,pgx)),(np-1)))); cnp=ep-sp+1;
	  plenv.xmin[1]=wl[sp]; plenv.xmax[1]=wl[ep];
	  plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
	  plenv.ymin[1]=-0.1*plenv.ymax[1];
	  plenv.ymax[2]=1.1*(MAX(0.5,(MIN(5.0,(fjmax(&(csq[sp]),cnp))))));
	  plenv.ymin[2]=-0.1*plenv.ymax[2];
	  plenv.ymax[3]=MAX(pgyo,pgy); plenv.ymin[3]=MIN(pgyo,pgy);
	  con_redo=1;
	}
	else con_redo=0; 
      }
      else if (!strncmp(pgch,"A",1)) {
	pgx=(MIN(wl[np-1],(MAX(pgx,wl[0]))));
	sp=MAX(0,(MIN((idxfval(wl,np,pgx)-cnp/2),(np-cnp)))); ep=sp+cnp-1;
	plenv.xmin[1]=wl[sp]; plenv.xmax[1]=wl[ep];
	plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
	plenv.ymin[1]=-0.1*plenv.ymax[1];
	plenv.ymax[2]=1.1*(MAX(0.5,(MIN(5.0,(fjmax(&(csq[sp]),cnp))))));
	plenv.ymin[2]=-0.1*plenv.ymax[2];
	plenv.ymax[3]=fjmax(&(fl[sp]),cnp);
	plenv.ymin[3]=MIN((-0.1*plenv.ymax[3]),(fjmin(&(fl[sp]),cnp)));
	con_redo=1;
      }
      else con_redo=0;
    }
    else if (pgy>cp->vpd3 && pgy<cp->vpu3 &&
	     (!strncmp(pgch,"A",1) ||
	      (!strncmp(pgch,"m",1) && plval && cact>-1 && (*act)[cact].act==COACT))) {
      cpgsvp(plenv.vpl,plenv.vpr,cp->vpd3,cp->vpu3);
      cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
      cpgsci(7); cpgslw(2.0*plenv.lw);
      cpgsch(0.9*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
      cpgsch(plenv.ch); cpgslw(plenv.lw); cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
      if (!strncmp(pgch,"?",1)) {
	fprintf(stderr,"\
Options for navigation/editing using detailed spectrum (flux panel):\n\
 Left mouse  : Re-centre panel at selected wavelength\n\
 Middle mouse: Re-centre panel (as above) and re-size window to suit\n\
                 new local extrema\n\
 Right mouse : Select custom window region.\n\
                 Right mouse again selects second corner of new window,\n\
                 any other entry to abort selection.\n\
 Space bar   : Write out to terminal information about cursor position.\n\
 b           : Define new bottom limit for plot\n\
 c           : Clip pixels from either combined or contributing spectra (act).\n\
                 c again selects second corner of clip window,\n\
               any other entry to abort selection.\n\
 l           : Define new left limit for plot\n\
 s           : Calculate normalized statistics over region (*)\n\
 r           : Define new right limit for plot\n\
 t           : Define new top limit for plot\n\
 u           : Unclip pixels from either combined or contributing spectra (act).\n\
                 u again selects second corner of clip window,\n\
                 any other entry to abort selection.\n\
 w           : Find equivalent width between this mark and next (*)\n\n");
	cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
      }

      if (!strncmp(pgch,"D",1)) {
	pgx=(MIN(wl[np-1],(MAX(pgx,wl[0]))));
	sp=MAX(0,(MIN((idxfval(wl,np,pgx)-cp->np/2),(np-cp->np))));
	ep=sp+cp->np-1; cnp=ep-sp+1;
	plenv.xmin[1]=wl[sp]; plenv.xmax[1]=wl[ep];
	plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
	plenv.ymin[1]=-0.1*plenv.ymax[1];
	plenv.ymax[2]=1.1*(MAX(0.5,(MIN(5.0,(fjmax(&(csq[sp]),cnp))))));
	plenv.ymin[2]=-0.1*plenv.ymax[2];
	plenv.ymax[3]=fjmax(&(fl[sp]),cnp);
	plenv.ymin[3]=MIN((-0.1*plenv.ymax[3]),(fjmin(&(fl[sp]),cnp)));
	con_redo=1;
      }
      else if (!strncmp(pgch,"X",1)) {
	strcpy(pgch,"\0"); cpgsci(2); pgxo=pgx; pgyo=pgy;
	cpgband(2,0,pgx,pgy,&pgx,&pgy,pgch);
	if (!strncmp(pgch,"X",1)) {
	  if (pgxo>pgx) { FSWAP(pgxo,pgx); }
	  if (pgyo>pgy) { FSWAP(pgyo,pgy); }
	  pgxo=(MAX(pgxo,wl[0])); pgx=(MIN(pgx,wl[np-1]));
	  sp=MAX(0,(MIN((idxfval(wl,np,pgxo)),(np-2))));
	  ep=MAX(0,(MIN((idxfval(wl,np,pgx)),(np-1)))); cnp=ep-sp+1;
	  plenv.xmin[1]=wl[sp]; plenv.xmax[1]=wl[ep];
	  plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
	  plenv.ymin[1]=-0.1*plenv.ymax[1];
	  plenv.ymax[2]=1.1*(MAX(0.5,(MIN(5.0,(fjmax(&(csq[sp]),cnp))))));
	  plenv.ymin[2]=-0.1*plenv.ymax[2];
	  plenv.ymax[3]=MAX(pgyo,pgy); plenv.ymin[3]=MIN(pgyo,pgy);
	  con_redo=1;
	}
      }
      else if (!strncmp(pgch,"A",1)) {
	pgx=(MIN(wl[np-1],(MAX(pgx,wl[0]))));
	sp=MAX(0,(MIN((idxfval(wl,np,pgx)-cnp/2),(np-cnp)))); ep=sp+cnp-1;
	plenv.xmin[1]=wl[sp]; plenv.xmax[1]=wl[ep];
	plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
	plenv.ymin[1]=-0.1*plenv.ymax[1];
	plenv.ymax[2]=1.1*(MAX(0.5,(MIN(5.0,(fjmax(&(csq[sp]),cnp))))));
	plenv.ymin[2]=-0.1*plenv.ymax[2];
	con_redo=1;
      }
      else if (!strncmp(pgch,"l",1)) {
	sp=MAX(0,(MIN((idxfval(wl,np,pgx)),(ep-1)))); cnp=ep-sp+1;
	plenv.xmin[1]=wl[sp];
	plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
	plenv.ymin[1]=-0.1*plenv.ymax[1];
	plenv.ymax[2]=1.1*(MAX(0.5,(MIN(5.0,(fjmax(&(csq[sp]),cnp))))));
	plenv.ymin[2]=-0.1*plenv.ymax[2];
	con_redo=1;
      }
      else if (!strncmp(pgch,"s",1)) {
	if (!plval || npcon==1) {
	  strcpy(pgch,"\0"); cpgsci(3); pgxo=pgx; pgyo=pgy;
	  cpgband(6,0,pgx,pgy,&pgx,&pgy,pgch);
	  if (!strncmp(pgch,"?",1)) {
	    fprintf(stderr,"Options for simple mathematics (flux panel):\n\
 m : Calculate statistics on normalized data between previous and this 'm'.\n\n");
	    cpgband(6,0,pgx,pgy,&pgx,&pgy,pgch);
	  }
	  if (!strncmp(pgch,"s",1)) {
	    if (pgxo>pgx) { FSWAP(pgxo,pgx); } if (pgyo>pgy) { FSWAP(pgyo,pgy); }
	    /* Find start and ending pixels on the screen */
	    csp=((csp=idxdval(&(cspec->rwl[sp]),cnp,pgxo))==-1) ? ep : csp+sp;
	    cep=((cep=idxdval(&(cspec->rwl[sp]),cnp,pgx))==-1) ? ep : cep+sp;
	    /* Decide whether calculating statistics is feasible given
	       current plot state */
	    if (plval && npcon==1) {
	      /** Calculate statistics for a single order **/
	      /* Find which spectrum and order we're dealing with */
	      i=0; while (!con[i][6]) i++;
	      /* Find start and ending pixels in relevant order */
	      cosp=(MAX(con[i][2],csp))-spec[con[i][0]].or[con[i][1]].csidx;
	      coep=(MIN(con[i][3],cep))-spec[con[i][0]].or[con[i][1]].csidx;
	      k=coep-cosp+1;
	      /* Allocate memory for data arrays */
	      if ((dat=darray(k))==NULL)
		errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor dat array of size %d",k);
	      if ((err=darray(k))==NULL)
		errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor err array of size %d",k);
	      if ((efl=darray(k))==NULL)
		errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor efl array of size %d",k);
	      if ((med=darray(k))==NULL)
		errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor med array of size %d",k);
	      /* Fill data arrays */
	      for (j=cosp,ndat=0; j<=coep; j++) {
		if (spec[con[i][0]].or[con[i][1]].rdst[j]==1 &&
		    spec[con[i][0]].or[con[i][1]].rdco[j]!=0.0) {
		  dat[ndat]=spec[con[i][0]].or[con[i][1]].rdfl[j]/
		    spec[con[i][0]].or[con[i][1]].rdco[j];
		  err[ndat]=spec[con[i][0]].or[con[i][1]].rder[j]/
		    spec[con[i][0]].or[con[i][1]].rdco[j];
		  efl[ndat]=spec[con[i][0]].or[con[i][1]].rdef[j]/
		    spec[con[i][0]].or[con[i][1]].rdco[j];
		  med[ndat++]=spec[con[i][0]].or[con[i][1]].rdme[j]/
		    spec[con[i][0]].or[con[i][1]].rdco[j];
		}
	      }
	      if (ndat) {
		/* Calculate statistics */
		if (!stats(dat,err,efl,med,NULL,ndat,2,&stdat))
		  errormsg("UVES_plot_cspec(): Error returned from stats()");
		pe=(ndat>1) ? gammq(0.5*(double)(ndat-1),0.5*stdat.chisq) : 0.0;
		/* Printf statistics to screen */
		fprintf(stderr,"Ntot=%6d Nval=%6d Spec=%d Order=%d\n",k,ndat,
			con[i][0]+1,con[i][1]+1);
		fprintf(stderr,"  m=%9.6lf +/- %8.6lf rms=%8.6lf med=%9.6lf\n",
			stdat.mean,stdat.emean,stdat.rms,stdat.med);
		fprintf(stderr,
		      " wm=%9.6lf +/- %8.6lf  me=%8.6lf chisqe=%8.6lf pe=%8.6lf\n",
			stdat.wmean,stdat.ewmean,stdat.meansig,stdat.rchisq,pe);
	      } else
		fprintf(stderr,"No valid pixels in this order in specified range");
	      /* Clean up */
	      free(dat); free(err); free(efl); free(med);
	    }
	    if (!plval) {	
	      /** Calculate statistics for a combined spectrum **/
	      k=cep-csp+1;
	      /* Find number of valid pixels */
	      for (j=csp,ndat=0; j<=cep; j++) if (cspec->st[j]==1) ndat++;
	      if (ndat) {
		/* Calculate statistics */
		if (!stats(&(cspec->no[csp]),&(cspec->ne[csp]),&(cspec->nf[csp]),
			   NULL,&(cspec->st[csp]),k,2,&stdat))
		  errormsg("UVES_plot_cspec(): Error returned from stats()");
		pe=(ndat>1) ? gammq(0.5*(double)(ndat-1),0.5*stdat.chisq) : 0.0;
		/* Printf statistics to screen */
		fprintf(stderr,"Ntot=%6d Nval=%6d\n",k,ndat);
		fprintf(stderr,"  m=%9.6lf +/- %8.6lf rms=%8.6lf med=%9.6lf\n",
			stdat.mean,stdat.emean,stdat.rms,stdat.med);
		fprintf(stderr,
		      " wm=%9.6lf +/- %8.6lf  me=%8.6lf chisqe=%8.6lf pe=%8.6lf\n",
			stdat.wmean,stdat.ewmean,stdat.meansig,stdat.rchisq,pe);
	      } else
		fprintf(stderr,
			"No valid pixels in combined spectrum in specified range");
	    }
	  }
	}
      }
      else if (!strncmp(pgch,"w",1)) {
	if (!plval || npcon==1) {
	  strcpy(pgch,"\0"); cpgsci(3); pgxo=pgx; pgyo=pgy;
	  cpgband(6,0,pgx,pgy,&pgx,&pgy,pgch);
	  if (!strncmp(pgch,"?",1)) {
	    fprintf(stderr,"Options for find equivalent width (flux panel):\n\
 w : Calculate equivalent width between previous and this 'w'.\n\n");
	    cpgband(6,0,pgx,pgy,&pgx,&pgy,pgch);
	  }
	  if (!strncmp(pgch,"w",1)) {
	    if (pgxo>pgx) { FSWAP(pgxo,pgx); } if (pgyo>pgy) { FSWAP(pgyo,pgy); }
	    /* Find start and ending pixels on the screen */
	    csp=((csp=idxdval(&(cspec->rwl[sp]),cnp,pgxo))==-1) ? ep : csp+sp;
	    cep=((cep=idxdval(&(cspec->rwl[sp]),cnp,pgx))==-1) ? ep : cep+sp;
	    /* Decide whether calculating EWs is feasible given
	       current plot state */
	    if (plval && npcon==1) {
	      /** Calculate EWs for a single order **/
	      /* Find which spectrum and order we're dealing with */
	      i=0; while (!con[i][6]) i++;
	      /* Find start and ending pixels in relevant order */
	      csp=(MAX(con[i][2],csp)); cep=(MIN(con[i][3],cep));
	      cosp=csp-spec[con[i][0]].or[con[i][1]].csidx;
	      coep=cep-spec[con[i][0]].or[con[i][1]].csidx; ndat=coep-cosp+1;
	      /* Allocate memory for data arrays */
	      if ((err=darray(ndat))==NULL)
		errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor err array of size %d",ndat);
	      /* Fill data arrays */
	      for (j=cosp,ndat=0,k=0; j<=coep; j++) {
		if (spec[con[i][0]].or[con[i][1]].rdst[j]==1 &&
		    spec[con[i][0]].or[con[i][1]].rdco[j]!=0.0) {
		  err[ndat++]=spec[con[i][0]].or[con[i][1]].rder[j]; k++;
		} else err[ndat++]=-INFIN;
	      }
	      /* Calculate EWs using the unwieghted method */
	      ew=eew=0.0;
	      cwl=0.5*(double)(pgx+pgxo); cwin=C_C_K*(double)(pgx-pgxo)/cwl;
	      if (!EW(&(cspec->wl[csp]),
		      &(spec[con[i][0]].or[con[i][1]].rdfl[cosp]),err,
		      &(spec[con[i][0]].or[con[i][1]].rdco[cosp]),ndat,&cwl,
		      cwin,par->disp,par->disp,NULL,NULL,NULL,&ew,&eew,NULL,NULL,
		      NULL,NULL,NBTOL,0,0))
		warnmsg("UVES_plot_cspec(): Error returned from EW() when\n\
\tfinding unweighted EW for section %f to %f of order %d in file\n\t%s",pgx,pgxo,
			con[i][1],spec[con[i][0]].file);
	      /* Print EWs to screen */
	      fprintf(stderr,"ws=%10.4f we=%10.4f Ntot=%6d Nval=%6d Spec=%d \
Order=%d\n",pgxo,pgx,ndat,k,con[i][0]+1,con[i][1]+1);
	      fprintf(stderr,"ew=%.6lg eew=%.6lg\n\n",ew,eew);
	      /* Clean up */
	      free(err);
	    }
	    if (!plval) {	
	      /** Calculate statistics for a combined spectrum **/
	      ndat=cep-csp+1;
	      /* Allocate memory for data arrays */
	      if ((err=darray(ndat))==NULL)
		errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor err array of size %d",ndat);
	      /* Fill data arrays */
	      for (j=csp,ndat=0,k=0; j<=cep; j++) {
		if (cspec->st[j]==1) { err[ndat++]=cspec->er[j]; k++; }
		else err[ndat++]=-INFIN;
	      }
	      /* Calculate EWs using the unwieghted method */
	      ew=eew=0.0;
	      cwl=0.5*(double)(pgx+pgxo); cwin=C_C_K*(double)(pgx-pgxo)/cwl;
	      if (!EW(&(cspec->wl[csp]),&(cspec->fl[csp]),err,&(cspec->co[csp]),
		      ndat,&cwl,cwin,par->disp,par->disp,NULL,NULL,NULL,&ew,&eew,
		      NULL,NULL,NULL,NULL,NBTOL,0,0))
		warnmsg("UVES_plot_cspec(): Error returned from EW() when\n\
\tfinding unweighted EW for section %f to %f in combined spectrum",pgx,pgxo);
	      /* Print EWs to screen */
	      fprintf(stderr,"ws=%10.4f we=%10.4f Ntot=%6d Nval=%6d\n",pgxo,pgx,
		      ndat,k);
	      fprintf(stderr,"ew=%.6lg eew=%.6lg\n\n",ew,eew);
	      /* Clean up */
	      free(err);
	    }
	  }
	}
      }
      else if (!strncmp(pgch,"r",1)) {
	ep=MIN((np-1),(MAX((idxfval(wl,np,pgx)),(sp+1)))); cnp=ep-sp+1;
	plenv.xmax[1]=wl[ep];
	plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
	plenv.ymin[1]=-0.1*plenv.ymax[1];
	plenv.ymax[2]=1.1*(MAX(0.5,(MIN(5.0,(fjmax(&(csq[sp]),cnp))))));
	plenv.ymin[2]=-0.1*plenv.ymax[2];
	con_redo=1;
      }
      else if (!strncmp(pgch,"t",1)) {
	plenv.ymax[3]=(pgy>plenv.ymin[3]) ? pgy : plenv.ymax[3];
	con_redo=0;
      }
      else if (!strncmp(pgch,"b",1)) {
	plenv.ymin[3]=(pgy<plenv.ymax[3]) ? pgy : plenv.ymin[3];
	con_redo=0;
      }
      else if (!strncmp(pgch," ",1)) {
	if (!plval || npcon==1) {
	  /* Find pixel on the screen */
	  csp=((csp=idxdval(&(cspec->rwl[sp]),cnp,pgx))==-1) ? ep : csp+sp;
	  if (plval && npcon==1) {
	    /** Write out cursor information for a single order **/
	    /* Find which spectrum and order we're dealing with */
	    i=0; while (!con[i][6]) i++;
	    /* Find pixel in relevant order */
	    cosp=(MAX(con[i][2],csp))-spec[con[i][0]].or[con[i][1]].csidx;
	    /* Write out relevant information for cursor position */
	    fprintf(stderr,"x=%12.6lf y=%-.8g\n",pgx,pgy);
	    fprintf(stderr,"w=%12.6lf Spec=%d Order=%d\n",cspec->wl[csp],
		    con[i][0]+1,con[i][1]+1);
	    fprintf(stderr,"Order   :  f=%-.8g  e=%-.8g  ef=%-.8g c=%9.6lf \
n=%6d\n",spec[con[i][0]].or[con[i][1]].rdfl[cosp],
		    spec[con[i][0]].or[con[i][1]].rder[cosp],
		    spec[con[i][0]].or[con[i][1]].rdef[cosp],
		    spec[con[i][0]].or[con[i][1]].rdco[cosp],cosp+1);
	    if (spec[con[i][0]].or[con[i][1]].rdco[cosp]!=0.0)
	      fprintf(stderr,"          nf=%9.6lf ne=%.8lg nef=%.8lg\n",
		      spec[con[i][0]].or[con[i][1]].rdfl[cosp]/
		      spec[con[i][0]].or[con[i][1]].rdco[cosp],
		      spec[con[i][0]].or[con[i][1]].rder[cosp]/
		      spec[con[i][0]].or[con[i][1]].rdco[cosp],
		      spec[con[i][0]].or[con[i][1]].rdef[cosp]/
		      spec[con[i][0]].or[con[i][1]].rdco[cosp]);
	    fprintf(stderr,"Combined:  f=%-.8g  e=%-.8g  ef=%-.8g c=%9.6lf \
n=%6d\n",cspec->fl[csp],cspec->er[csp],cspec->ef[csp],cspec->co[csp],csp+1);
	    fprintf(stderr,"          nf=%9.6lf ne=%.8lg nef=%.8lg\n",
		    cspec->no[csp],cspec->ne[csp],cspec->nf[csp]);
	    if (cspec->st[csp]==1) {
	      fprintf(stderr," N(bef)=%d chisq(bef)=%9.6lf p(bef)=%8.6lf\n",
		      cspec->ncb[csp],cspec->csq[csp],
		      gammq(0.5*(double)(cspec->ncb[csp]-1),0.5*cspec->csq[csp]*
			    (double)(cspec->ncb[csp]-1)));
	      fprintf(stderr," N(aft)=%d chisq(aft)=%9.6lf p(aft)=%8.6lf\n\n",
		      cspec->nccb[csp],cspec->ccsq[csp],
		      gammq(0.5*(double)(cspec->nccb[csp]-1),0.5*cspec->ccsq[csp]*
			    (double)(cspec->nccb[csp]-1)));
	    }
	  }	    
	  if (!plval) {	
	    /** Write out cursor information for combined spectrum **/
	    fprintf(stderr,"x=%12.6lf y=%-.8g\n",pgx,pgy);
	    fprintf(stderr,"w=%12.6lf f=%-.8g  e=%-.8g  ef=%-.8g c=%9.6lf n=%6d\n",
		    cspec->wl[csp],cspec->fl[csp],cspec->er[csp],cspec->ef[csp],
		    cspec->co[csp],csp+1);
	    fprintf(stderr,"              nf=%9.6lf ne=%.8lg nef=%.8lg\n",
		    cspec->no[csp],cspec->ne[csp],cspec->nf[csp]);
	    if (cspec->st[csp]==1) {
	      fprintf(stderr,"N(bef)=%d chisq(bef)=%9.6lf p(bef)=%8.6lf\n",
		      cspec->ncb[csp],cspec->csq[csp],
		      gammq(0.5*(double)(cspec->ncb[csp]-1),0.5*cspec->csq[csp]*
			    (double)(cspec->ncb[csp]-1)));
	      fprintf(stderr,"N(aft)=%d chisq(aft)=%9.6lf p(aft)=%8.6lf\n\n",
		      cspec->nccb[csp],cspec->ccsq[csp],
		      gammq(0.5*(double)(cspec->nccb[csp]-1),0.5*cspec->ccsq[csp]*
			    (double)(cspec->nccb[csp]-1)));
	    }
	  }
	} else fprintf(stderr,"x=%12.6lf y=%-.8g\n\n",pgx,pgy);
	con_redo=0;
      }
      else if (!strncmp(pgch,"c",1) ||
	       (!strncmp(pgch,"m",1) && plval && cact>-1 && (*act)[cact].act==COACT)) {
	/* Clip pixels */
	cpgsci(3);
	if (strncmp(pgch,"m",1)) {
	  strcpy(pgch,"\0"); pgxo=pgx; pgyo=pgy; cpgband(2,0,pgx,pgy,&pgx,&pgy,pgch);
	} else {
	  pgxo=(*act)[cact].d[0]; pgyo=(*act)[cact].d[2];
	  pgx=(*act)[cact].d[1]; pgy=(*act)[cact].d[3];
	  strcpy(pgch,"c\0");
	}
	if (!strncmp(pgch,"c",1)) {
	  if (pgxo>pgx) { FSWAP(pgxo,pgx); } if (pgyo>pgy) { FSWAP(pgyo,pgy); }
	  csp=((csp=idxdval(&(cspec->rwl[sp]),cnp,pgxo))==-1) ? ep : csp+sp;
	  cep=((cep=idxdval(&(cspec->rwl[sp]),cnp,pgx))==-1) ? ep : cep+sp;
	  if (plval) {
	    /* Clip pixels from individual orders */
	    for (i=0,k=0; i<ncon; i++) {
	      if (con[i][6] && (con[i][2]<=cep && con[i][3]>=csp)) {
		cosp=(MAX(con[i][2],csp))-spec[con[i][0]].or[con[i][1]].csidx;
		coep=(MIN(con[i][3],cep))-spec[con[i][0]].or[con[i][1]].csidx;
		for (j=cosp,selpix=0; j<=coep; j++) {
		  if (spec[con[i][0]].or[con[i][1]].rdfl[j]<pgy &&
		      spec[con[i][0]].or[con[i][1]].rdfl[j]>pgyo) {
		    spec[con[i][0]].or[con[i][1]].ordst[j]=
		      spec[con[i][0]].or[con[i][1]].rdst[j];
		    spec[con[i][0]].or[con[i][1]].rdst[j]=CCLIP; selpix++;
		  }
		}
		/* Remember action by filling in an action report */
		if (selpix) {
		  if (!(*nact)) {
		    if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
		      errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
		  }
		  else
		    if (!(*act=(action *)
			  realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
		      errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
		  /* If in paused replay mode then shift later actions to end */
		  if (par->replay && rp->exit==2) {
		    for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
		    cact=rp->cact++; (*nact)++;
		  } else cact=(*nact)++;
		  (*act)[cact].act=COACT; (*act)[cact].nxyp=0;
		  (*act)[cact].d[0]=pgxo; (*act)[cact].d[1]=pgx;
		  (*act)[cact].d[2]=pgyo; (*act)[cact].d[3]=pgy;
		  (*act)[(cact)].i[0]=con[i][0]; (*act)[(cact)].i[1]=con[i][1];
		  (*act)[cact].rcmb=0; (*act)[cact].val=1; k++;
		}
	      }
	    }
	    /* If some new actions were added, record how many */
	    if (*nact) (*act)[cact].nordact=k;
	  }
	  else {
	    /* Clip pixels from combined spectrum */
	    for (j=csp,selpix=0; j<=cep; j++) {
	      if (cspec->fl[j]<pgy && cspec->fl[j]>pgyo) {
		cspec->ost[j]=cspec->st[j]; cspec->st[j]=CCLIP; selpix++;
		cspec->no[j]=1.0; cspec->ne[j]=cspec->nf[j]=-INFIN;
		no[j]=1.0; ne[j]=-1.0;
	      }
	    }
	    /* Remember action by filling in an action report */
	    if (selpix) {
	      if (!(*nact)) {
		if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
	      }
	      else
		if (!(*act=(action *)
		      realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
	      /* If in paused replay mode then shift later actions to end */
	      if (par->replay && rp->exit==2) {
		for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
		cact=rp->cact++; (*nact)++;
	      } else cact=(*nact)++;
	      (*act)[cact].act=CCACT; (*act)[cact].nxyp=0; (*act)[cact].d[0]=pgxo;
	      (*act)[cact].d[1]=pgx; (*act)[cact].d[2]=pgyo; (*act)[cact].d[3]=pgy;
	      (*act)[cact].rcmb=0; (*act)[cact].nordact=1; (*act)[cact].val=1; 
	    }
	  }
	}
	con_redo=0; cp->rccsp=1;
      }
      else if (!strncmp(pgch,"u",1)) {
	/* Unclip pixels */
	strcpy(pgch,"\0"); cpgsci(3); pgxo=pgx; pgyo=pgy;
	cpgband(2,0,pgx,pgy,&pgx,&pgy,pgch);
	if (!strncmp(pgch,"u",1)) {
	  if (pgxo>pgx) { FSWAP(pgxo,pgx); } if (pgyo>pgy) { FSWAP(pgyo,pgy); }
	  csp=((csp=idxdval(&(cspec->rwl[sp]),cnp,pgxo))==-1) ? ep : csp+sp;
	  cep=((cep=idxdval(&(cspec->rwl[sp]),cnp,pgx))==-1) ? ep : cep+sp;
	  if (plval) {
	    /* Unclip pixels from individual orders */
	    for (i=0,k=0; i<ncon; i++) {
	      if (con[i][6] && (con[i][2]<=cep && con[i][3]>=csp)) {
		cosp=(MAX(con[i][2],csp))-spec[con[i][0]].or[con[i][1]].csidx;
		coep=(MIN(con[i][3],cep))-spec[con[i][0]].or[con[i][1]].csidx;
		for (j=cosp,selpix=0; j<=coep; j++) {
		  if (spec[con[i][0]].or[con[i][1]].rdfl[j]<pgy &&
		      spec[con[i][0]].or[con[i][1]].rdfl[j]>pgyo) {
		    spec[con[i][0]].or[con[i][1]].ordst[j]=
		      spec[con[i][0]].or[con[i][1]].rdst[j];
		    spec[con[i][0]].or[con[i][1]].rdst[j]=1; selpix++;
		  }
		}
		/* Remember action by filling in an action report */
		if (selpix) {
		  if (!(*nact)) {
		    if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
		      errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
		  }
		  else
		    if (!(*act=(action *)
			  realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
		      errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
		  /* If in paused replay mode then shift later actions to end */
		  if (par->replay && rp->exit==2) {
		    for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
		    cact=rp->cact++; (*nact)++;
		  } else cact=(*nact)++;
		  (*act)[cact].act=UOACT; (*act)[cact].nxyp=0;
		  (*act)[cact].d[0]=pgxo; (*act)[cact].d[1]=pgx;
		  (*act)[cact].d[2]=pgyo; (*act)[cact].d[3]=pgy;
		  (*act)[(cact)].i[0]=con[i][0]; (*act)[(cact)].i[1]=con[i][1];
		  (*act)[cact].rcmb=0; (*act)[cact].val=1; k++;
		}
	      }
	    }
	    /* If some new actions were added, record how many */
	    if (*nact) (*act)[cact].nordact=k;
	  }
	  else {
	    /* Unclip pixels from combined spectrum */
	    for (j=csp,selpix=0; j<=cep; j++) {
	      if (cspec->fl[j]<pgy && cspec->fl[j]>pgyo) {
		cspec->ost[j]=cspec->st[j]; cspec->st[j]=1; selpix++;
		cspec->no[j]=cspec->fl[j]/cspec->co[j];
		cspec->ne[j]=cspec->er[j]/cspec->co[j];
		cspec->nf[j]=cspec->ef[j]/cspec->co[j];
		no[j]=(float)cspec->no[j]; ne[j]=(float)cspec->ne[j];
	      }
	    }
	    /* Remember action by filling in an action report */
	    k=0; if (selpix) {
	      if (!(*nact)) {
		if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
	      }
	      else
		if (!(*act=(action *)
		      realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
	      /* If in paused replay mode then shift later actions to end */
	      if (par->replay && rp->exit==2) {
		for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
		cact=rp->cact++; (*nact)++;
	      } else cact=(*nact)++;
	      (*act)[cact].act=UCACT; (*act)[cact].nxyp=0; (*act)[cact].d[0]=pgxo;
	      (*act)[cact].d[1]=pgx; (*act)[cact].d[2]=pgyo; (*act)[cact].d[3]=pgy;
	      (*act)[cact].rcmb=0; (*act)[cact].val=1; k++;
	    }
	    /* If some new actions were added, record how many */
	    if (*nact) (*act)[cact].nordact=k;
	  }
	}
	con_redo=0; cp->rccsp=1;
      }
      else con_redo=0;
    }
    else if (!strncmp(pgch,"y",1)) {
      plenv.ymax[3]=1.1*fjmax(&(fl[sp]),cnp);
      plenv.ymin[3]=MIN((-0.1*plenv.ymax[3]),(fjmin(&(fl[sp]),cnp)));
      con_redo=0;
    }
    else if (!strncmp(pgch,"s",1)) {
      if (!UVES_wUPLfile(cspec->UPLfile,spec,nspec,*act,*nact,cspec,par))
	errormsg("Unknown error returned from UVES_wUPLfile()");
      *nact_save=*nact;
      if (!UVES_wFITSfile(spec,nspec,cspec,par))
	errormsg("Unknown error returned from UVES_wFITSfile()");
      if (par->dat && !UVES_wDATfile(cspec->DATfile,cspec))
	errormsg("Uknown error returned from UVES_wDATfile()");
      con_redo=0;
    }
    else if (!strncmp(pgch,".",1)) {
      ftemp=0.5*(plenv.ymax[3]-plenv.ymin[3]);
      plenv.ymax[3]+=ftemp; plenv.ymin[3]-=ftemp;
      con_redo=0;
    }
    else if (!strncmp(pgch,",",1)) {
      ftemp=0.25*(plenv.ymax[3]-plenv.ymin[3]);
      plenv.ymax[3]-=ftemp; plenv.ymin[3]+=ftemp;
      con_redo=0;
    }
    else if (!strncmp(pgch,">",1)) {
      sp=MAX(0,sp-cnp/2); ep=MIN(np-1,ep+cnp/2); cnp=ep-sp+1; 
      plenv.xmin[1]=wl[sp]; plenv.xmax[1]=wl[ep];
      plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
      plenv.ymin[1]=-0.1*plenv.ymax[1];
      con_redo=1;
    }
    else if (!strncmp(pgch,"<",1)) {
      sp+=cnp/4; ep-=cnp/4; cnp=ep-sp+1; 
      plenv.xmin[1]=wl[sp]; plenv.xmax[1]=wl[ep];
      plenv.ymax[1]=1.1*(MAX(1.0,(fjmax(&(ncb[sp]),cnp))));
      plenv.ymin[1]=-0.1*plenv.ymax[1];
      con_redo=1;
    }
    else if (!strncmp(pgch,"k",1)) {
      if (!plval || tscon || !ncon) tscon=nscon=0;
      else {
	tscon=1;
	for (i=0,nscon=0; i<ncon; i++)
	  if (con[i][6]) { scon[nscon][0]=con[i][0]; scon[nscon++][1]=con[i][1]; }
      }
      con_redo=1;
    }
    else if (!strncmp(pgch,"o",1) && !par->combmeth) {
      if (!par->scalmeth) cp->irank=(PMOD(cp->irank,1,1,0));
      else if (par->scalmeth==1) cp->irank=(PMOD(cp->irank,1,2,0));
      con_redo=1;
    }
    else if (!strncmp(pgch,"v",1)) {
      plval=(plval) ? 0 : 1; con_redo=1;
    }
    else if (plval &&
	     (pgy>plenv.vpu || !strncmp(pgch,"n",1) || !strncmp(pgch,"p",1))) {
      /* Select which order(s) to plot */
      step=0.9*(plenv.vpr-plenv.vpl)/(double)ncon;
      idx=(pgx<plenv.vpl+0.1*(plenv.vpr-plenv.vpl)) ? -1 :
	(int)((pgx-plenv.vpl-0.1*(plenv.vpr-plenv.vpl))/step);
      if (idx>=0) {
	if (idx>=ncon) idx=ncon-1;
	if (!strncmp(pgch,"n",1) && npcon) {
	  j=idx; while (!con[j][6]) j=(PMOD(j,1,(ncon-1),0));
	  con[j][6]=0; j=(PMOD(j,1,(ncon-1),0)); con[j][6]=1;
	  con_redo=0;
	}
	else if (!strncmp(pgch,"p",1) && npcon) {
	  j=idx; while (!con[j][6]) j=(j) ? j-1 : ncon-1;
	  con[j][6]=0; j=(j) ? j-1 : ncon-1; con[j][6]=1;
	  con_redo=0;
	}
	else if (!strncmp(pgch,"A",1)) {
	  for (i=0; i<ncon; i++) con[i][6]=0; con[idx][6]=1;
	  con_redo=0;
	}
	else if (!strncmp(pgch,"D",1)) {
	  for (i=0; i<idx; i++) con[i][6]=0;
	  for (i=idx; i<ncon; i++) con[i][6]=1;
	  con_redo=0;
	}
	else if (!strncmp(pgch,"X",1)) {
	  con[idx][6]=(con[idx][6]) ? 0 : 1;
	  con_redo=0;
	}
	else if (!strncmp(pgch,"i",1)) {
	  fprintf(stderr,"Spec=%d Order=%d\n\n",con[idx][0]+1,con[idx][1]+1);
	}
	else if (!strncmp(pgch,"l",1)) {
	  for (i=0; i<=idx; i++) con[i][6]=1;
	  for (i=idx+1; i<ncon; i++) con[i][6]=0;
	  con_redo=0;
	}
	else if (!strncmp(pgch,"e",1)) {
	  tscon=1;
	  for (i=0,nscon=0; i<ncon; i++) {
	    if (con[i][6] || i==idx) {
	      scon[nscon][0]=con[i][0]; scon[nscon++][1]=con[i][1];
	    }
	  }
	  for (j=0; j<spec[con[idx][0]].nor; j++) {
	    if (spec[con[idx][0]].or[j].nuse>=MINUSE && j!=con[idx][1]) {
	      scon[nscon][0]=con[idx][0]; scon[nscon++][1]=j;
	    }
	  }
	  con_redo=1;
	}
	else if (!strncmp(pgch,"x",1)) {
	  tscon=1;
	  for (i=0,nscon=0; i<nspec; i++) {
	    if (i!=con[idx][0]) {
	      for (j=0; j<spec[i].nor; j++) {
		scon[nscon][0]=i; scon[nscon++][1]=j;
	      }
	    }   
	  }
	  con_redo=1;
	}
	else if (!strncmp(pgch,"2",1) || !strncmp(pgch,"1",1) ||
		 !strncmp(pgch,"3",1) || !strncmp(pgch,"4",1) ||
		 (!strncmp(pgch,"m",1) && cact>-1 && (*act)[cact].act==SOACT &&
		  (con[idx][0]!=(*act)[cact].i[0] || con[idx][1]!=(*act)[cact].i[1])))
	  {
	  /* Manually scale a particular order */
	  if (!strncmp(pgch,"3",1)) scale=1.0+sscale;
	  else if (!strncmp(pgch,"4",1)) scale=1.0+bscale;
	  else if (!strncmp(pgch,"2",1)) scale=1.0-sscale;
	  else if (!strncmp(pgch,"1",1)) scale=1.0-bscale;
	  else scale=(*act)[cact].d[0];
	  for (i=0; i<spec[con[idx][0]].or[con[idx][1]].nrdp; i++) {
	    spec[con[idx][0]].or[con[idx][1]].rdfl[i]*=scale;
	    spec[con[idx][0]].or[con[idx][1]].rder[i]*=scale;
	    spec[con[idx][0]].or[con[idx][1]].rdef[i]*=scale;
	    spec[con[idx][0]].or[con[idx][1]].rdme[i]*=scale;
	  }
	  /* If the same spectrum/order and action as last action, simply
	     update last action with new scale factor */
	  if (*nact) {
	    k=0; if ((*act)[cact].act==SOACT && (*act)[cact].i[0]==con[idx][0] &&
		     (*act)[cact].i[1]==con[idx][1]) k=1;
	    if (par->replay && rp->exit==2) { if (sclinit!=cact) k=0; }
	    if (k) (*act)[cact].d[0]*=scale;
	    else {
	      if (!(*act=(action *)
		    realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
		errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
	      /* If in paused replay mode then shift later actions to end */
	      if (par->replay && rp->exit==2) {
		for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
		cact=rp->cact++; (*nact)++; sclinit=cact;
	      } else cact=(*nact)++;
	      (*act)[cact].act=SOACT; (*act)[cact].nxyp=0; (*act)[cact].d[0]=scale;
	      (*act)[cact].i[0]=con[idx][0]; (*act)[cact].i[1]=con[idx][1];
	      (*act)[cact].rcmb=0; (*act)[cact].nordact=(*act)[cact].val=1;
	    }
	  }
	  else {
	    if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
	      errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
	    cact=0;
	    (*act)[cact].act=SOACT; (*act)[cact].nxyp=0; (*act)[cact].d[0]=scale;
	    (*act)[cact].i[0]=con[idx][0]; (*act)[cact].i[1]=con[idx][1];
	    (*act)[cact].rcmb=0; (*act)[cact].nordact=(*act)[cact].val=1;
	    (*nact)++;
	  }
	  /* Record scaling for the relevant order's undo information */
	  if (*nact) {
	    spec[con[idx][0]].or[con[idx][1]].oscl=
	      spec[con[idx][0]].or[con[idx][1]].scl;
	    spec[con[idx][0]].or[con[idx][1]].scl*=scale;
	  }
	  con_redo=0; cp->rccsp=cp->rscsp=1;
	}
      }
    }
    else if (!strncmp(pgch,"i",1)) {
      /*** Interpolate a new continuum, either to the final combined
	  spectrum or to a single order (if only one plotted) ***/
      /* Decide whether interpolation is feasible given current plot state */
      if (!plval || npcon==1) {
	/* Refresh number of pairs of xy pairs of positions */
	if (plval && npcon==1) {
	  /** Fit new continuum for a single order **/
	  /* Find which spectrum and order we're dealing with */
	  i=0; while (!con[i][6]) i++;
	  /* Find initial fit region */
	  if (par->replay && rp->exit==3) {
	    /* If in replay mode we must use the values stores in the action */
	    if ((fsp=idxdval(&(cspec->rwl[sp]),cnp,(*act)[rp->cact].d[0]))==-1) {
	      nferrormsg("UVES_plot_cspec(): Cannot find pixel with\n\
\twavelength %lf in combined spectrum",(*act)[rp->cact].d[0]); return 0;
	    } else fsp+=sp;
	    if ((fep=idxdval(&(cspec->rwl[fsp]),cnp+sp-fsp,(*act)[rp->cact].d[1]))
		==-1) fep=ep;
	    else fep+=fsp;
	    if ((fbsp=idxdval(&(cspec->rwl[fsp]),cnp+sp-fsp,(*act)[rp->cact].d[2]))
		==-1) fbsp=fep-2;
	    else fbsp+=fsp;
	    if ((fbep=idxdval(&(cspec->rwl[fbsp]),cnp+sp-fbsp,
			      (*act)[rp->cact].d[3]))==-1) fbep=fep-1;
	    else fbep+=fbsp;
	    fnp=fep-fsp+1;
	    if (fnp<4) {
	      nferrormsg("UVES_plot_cspec(): Too few (=%d) pixels\n\
\tto interpolate across"); return 0;
	    }
	    fbsp=(MAX(fbsp,(fsp+1))); fbep=(MAX(fbep,(fbsp+1)));
	    fep=(MAX(fep,(fbep+1))); fbsnp=fbsp-fsp+1; fbenp=fep-fbep+1;
	    /* Initialize number of x-y pairs */
	    nxyp=(*act)[rp->cact].nxyp;
	  } else {
	    fsp=con[i][2]+(int)(FFRAC*con[i][4]);
	    fep=con[i][3]-(int)(FFRAC*con[i][4]); fnp=fep-fsp+1;
	    fbsp=fsp+(int)(FBFRAC*(double)fnp); fbsnp=fbsp-fsp+1;
	    fbep=fep-(int)(FBFRAC*(double)fnp); fbenp=fep-fbep+1;
	    nxyp=2;
	  }
	  fosp=fsp-spec[con[i][0]].or[con[i][1]].csidx; foep=fosp+fnp-1;
          /* Allocate memory for fitting arrays */
          if ((fco=darray(con[i][4]))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fco array of size %d",con[i][4]);
          if ((fx=darray(nxyp))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fx array of size %d",nxyp);
          if ((fy=darray(nxyp))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fy array of size %d",nxyp);
          if ((fy2=darray(nxyp))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fy2 array of size %d",nxyp);
	  if (par->replay && rp->exit==3) {
	    /* Copy over node information and free x-y pair array in
	       original action */
	    for (i=0; i<nxyp; i++) {
	      fx[i]=(*act)[rp->cact].xyp[i].x1; fy[i]=(*act)[rp->cact].xyp[i].y1;
	    }
	    free((*act)[rp->cact].xyp); (*act)[rp->cact].xyp=NULL;
	  } else {
	    /* Initialize y end points of node array to reasonable values */
	    fy[0]=spec[con[i][0]].or[con[i][1]].rdco[fosp];
	    fy[1]=spec[con[i][0]].or[con[i][1]].rdco[foep];
	  }
	  /* Fill in fitted continuum array */
	  for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx,l=con[i][2];
	       j<con[i][4]; j++,k++,l++)
	    fco[j]=spec[con[i][0]].or[con[i][1]].rdco[k];
	  /* Start fit selection loop */
	  strcpy(pgch,"\0");
	  while (strncmp(pgch,"q",1) && strncmp(pgch,"a",1)) {
	    /* End points of spline node array are always ends of fitting region */
	    fx[0]=cspec->wl[fsp]; fx[nxyp-1]=cspec->wl[fep];
	    fosp=fsp-spec[con[i][0]].or[con[i][1]].csidx; foep=fosp+fnp-1;
	    /* Plot raw data again */
	    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd3,cp->vpu3); cpgsci(0);
	    cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
	    cpgsci(15);
	    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
	    cpgsci(4); cpgslw(2.0*plenv.lw);
	    cpgsch(0.9*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
	    cpgsch(0.8*plenv.ch); cpgslw(plenv.lw); cpgsci(3);
	    cpglab(" ",plenv.ylab[3]," "); cpgsci(con[i][5]);
	    for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
		 j<con[i][4]; j++,k++)
	      or[j]=spec[con[i][0]].or[con[i][1]].rdfl[k];
	    cpgbin(con[i][4],&(wl[con[i][2]]),or,1);
	    cpgsci(5); cpgsls(4);
	    cpgmove(plenv.xmin[1],0.0); cpgdraw(plenv.xmax[1],0.0); cpgsls(1);
	    /* Reinitalize fitted continuum array */
	    for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx,l=con[i][2];
	       j<con[i][4]; j++,k++,l++)
	    fco[j]=spec[con[i][0]].or[con[i][1]].rdco[k];
	    /* Calculate spline */
	    if (!spline(fx,fy,nxyp,0.0,0.0,fy2,1))
	      errormsg("UVES_plot_cspec(): Error returned from spline()");
	    for (j=fsp,k=fsp-con[i][2]; j<=fep; j++,k++)
	      fco[k]=splint(fx,fy,fy2,nxyp,cspec->wl[j]);
	    /* Blend two continua together at ends over the blend region */
	    for (j=0,k=fsp-con[i][2],l=fosp; j<fbsnp; j++,k++,l++) {
	      frac=(double)(j+1)/(double)fbsnp;
	      fco[k]=frac*fco[k]+(1.0-frac)*spec[con[i][0]].or[con[i][1]].rdco[l];
	    }
	    for (j=0,k=fep-con[i][2],l=foep; j<fbenp; j++,k--,l--) {
	      frac=(double)(j+1)/(double)fbenp;
	      fco[k]=frac*fco[k]+(1.0-frac)*spec[con[i][0]].or[con[i][1]].rdco[l];
	    }
	    /* Plot the original continuum over the data */
	    cpgsci((PMOD(con[i][5],1,14,2))); cpgline(cnp,&(wl[sp]),&(co[sp]));
	    /* Plot new continuum */
	    for (j=0,k=fsp-con[i][2]; j<fnp; j++,k++) or[j]=fco[k];
	    cpgsci((PMOD(con[i][5],2,14,2))); cpgline(fnp,&(wl[fsp]),or);
	    /* Mark spline nodes */
	    cpgsch(0.6*plenv.ch);
	    for (j=0; j<nxyp; j++) cpgpt1((float)fx[j],(float)fy[j],8);
	    /* Plot fit and blend region bounds */
	    cpgsci(3); cpgsls(2);
	    cpgmove(wl[fbsp],plenv.ymin[3]); cpgdraw(wl[fbsp],plenv.ymax[3]);
	    cpgmove(wl[fbep],plenv.ymin[3]); cpgdraw(wl[fbep],plenv.ymax[3]);
	    cpgsls(1);
	    cpgmove(wl[fsp],plenv.ymin[3]); cpgdraw(wl[fsp],plenv.ymax[3]);
	    cpgmove(wl[fep],plenv.ymin[3]); cpgdraw(wl[fep],plenv.ymax[3]);
	    /* Alter fit by selecting new bounds or selecting/deselecting
	       spline nodes points */
	    cpgsci((PMOD(con[i][5],3,14,2)));
	    cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	    if (!strncmp(pgch,"?",1)) {
	      fprintf(stderr,"\
Options for windowing new plotting region:\n\
 Left mouse  : Select left, right, or smoothing limit to change.\n\
                 Any key to place limit at new position and 'r'\n\
                 registers selected smoothing limit to corresponding\n\
                 fitting limit.\n\
 Right mouse : Remove spline node nearest cursor\n\
 Middle mouse: Move spline node nearest cursor.\n\
                 Any key to place at new position.\n\
 Space       : New spline node at cursor\n\
 a           : Accept fit\n\
 q           : Quit fitting procedure\n\n");
	      cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	    }
	    if (!strncmp(pgch,"A",1)) {
	      clickx=0.01*(plenv.xmax[1]-plenv.xmin[1]);
	      clicknp=idxfval(&(wl[con[i][2]]),con[i][4],wl[con[i][2]]+clickx);
	      clicknp+=2;
	      if (pgx>wl[fsp]-clickx && pgx<wl[fsp]+clickx) {
		cpgsci(15);
		cpgmove(wl[fsp],plenv.ymin[3]); cpgdraw(wl[fsp],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		k=idxfval(&(wl[fsp]),fnp,(float)fx[1])-1; k+=fsp;
		if (!strncmp(pgch,"r",1)) {
		  fsp=(MIN(k,(fbsp-1))); fnp=fep-fsp+1; fbsnp=1;
		}
		else {
		  fsp=((fsp=idxfval(&(wl[con[i][2]]),con[i][4],pgx))==-1) ?
		    fbep-clicknp : fsp+con[i][2]; fsp=MIN(fsp,(fbep-clicknp));
		  fsp=(MIN(k,fsp)); fnp=fep-fsp+1;
		  cpgsci(15); cpgmove(wl[fbsp],plenv.ymin[3]);
		  cpgdraw(wl[fbsp],plenv.ymax[3]);
		  fbsp=MAX(fbsp,(MIN((fbep-1),(fsp+fbsnp-1)))); fbsnp=fbsp-fsp+1;
		}
		fosp=fsp-spec[con[i][0]].or[con[i][1]].csidx; foep=fosp+fnp-1;
		fy[0]=spec[con[i][0]].or[con[i][1]].rdco[fosp];
	      }
	      else if (pgx>wl[fep]-clickx && pgx<wl[fep]+clickx) {
		cpgsci(15);
		cpgmove(wl[fep],plenv.ymin[3]); cpgdraw(wl[fep],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		k=idxfval(&(wl[fsp]),fnp,(float)fx[nxyp-2]); k+=fsp;
		if (!strncmp(pgch,"r",1)) {
		  fep=(MAX(k,(fbep+1))); fnp=fep-fsp+1; fbenp=1;
		}
		else {
		  fep=((fep=idxfval(&(wl[con[i][2]]),con[i][4],pgx))==-1) ?
		    con[i][3] : fep+con[i][2]; fep=MAX(fep,(fbsp+clicknp));
		  fep=(MAX(k,fep)); fnp=fep-fsp+1;
		  cpgsci(15); cpgmove(wl[fbep],plenv.ymin[3]);
		  cpgdraw(wl[fbep],plenv.ymax[3]);
		  fbep=MIN(fbep,(MAX((fbsp+1),(fep-fbenp+1)))); fbenp=fep-fbep+1;
		}
		fosp=fsp-spec[con[i][0]].or[con[i][1]].csidx; foep=fosp+fnp-1;
		fy[nxyp-1]=spec[con[i][0]].or[con[i][1]].rdco[foep];
	      }
	      else if (pgx>wl[fbsp]-clickx && pgx<wl[fbsp]+clickx) {
		cpgsci(15); cpgmove(wl[fbsp],plenv.ymin[3]);
		cpgdraw(wl[fbsp],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fbsp=fsp+1; fbsnp=1; }
		else {
		  fbsp=((fbsp=idxfval(&(wl[con[i][2]]),con[i][4],pgx))==-1) ?
		    fbep-1 : fbsp+con[i][2];
		  fbsp=MAX((fsp+1),(MIN(fbsp,(fbep-1)))); fbsnp=fbsp-fsp+1;
		}
	      }
	      else if (pgx>wl[fbep]-clickx && pgx<wl[fbep]+clickx) {
		cpgsci(15); cpgmove(wl[fbep],plenv.ymin[3]);
		cpgdraw(wl[fbep],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fbep=fep-1; fbenp=1; }
		else {
		  fbep=((fbep=idxfval(&(wl[con[i][2]]),con[i][4],pgx))==-1) ?
		    fep-1 : fbep+con[i][2];
		  fbep=MIN((fep-1),(MAX(fbep,(fbsp+1)))); fbenp=fep-fbep+1;
		}
	      }
	    }
	    else if (!strncmp(pgch,"X",1) && nxyp>2) {
	      /* Remove spline node nearest cursor */
	      k=1; cdist=(pgx-(float)fx[k])/(plenv.xmax[1]-plenv.xmin[1]);
	      cdist=sqrt(cdist*cdist);
	      for (j=2; j<nxyp-1; j++) {
		ftemp=(pgx-(float)fx[j])/(plenv.xmax[1]-plenv.xmin[1]);
		ftemp=sqrt(ftemp*ftemp); if (ftemp<cdist) { cdist=ftemp; k=j; }
	      }
	      nxyp--; for (j=k; j<nxyp; j++) { fx[j]=fx[j+1]; fy[j]=fy[j+1]; }
	      /* Reallocate memory in accordance with removal */
	      if (!(fx=(double *)realloc(fx,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fx\n\
\tarray of size %d",nxyp);
	      if (!(fy=(double *)realloc(fy,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fy\n\
\tarray of size %d",nxyp);
	      if (!(fy2=(double *)realloc(fy2,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fy2\n\
\tarray of size %d",nxyp);
	    }
	    else if (!strncmp(pgch,"D",1)) {
	      /* Move spline node nearest cursor */
	      k=0; cdist=(pgx-(float)fx[k])/(plenv.xmax[1]-plenv.xmin[1]);
	      cdist=sqrt(cdist*cdist);
	      for (j=1; j<nxyp; j++) {
		ftemp=(pgx-(float)fx[j])/(plenv.xmax[1]-plenv.xmin[1]);
		ftemp=sqrt(ftemp*ftemp); if (ftemp<cdist) { cdist=ftemp; k=j; }
	      }
	      /* Erase symbol for selected node */
	      cpgsci(15); cpgsch(0.6*plenv.ch);
	      cpgpt1((float)fx[k],(float)fy[k],8);
	      cpgsci((PMOD(con[i][5],4,14,2)));
	      cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	      if (k!=0 && k!=nxyp-1)
		fx[k]=(MAX((MIN(((double)pgx),((1.0-1.e-7)*fx[k+1]))),
			   ((1.0+1.e-7)*fx[k-1])));
	      fy[k]=(double)pgy;
	    }
	    else if (!strncmp(pgch," ",1) && pgx>(float)fx[0] &&
		     pgx<(float)fx[nxyp-1]) {
	      /* Find index at which to insert new node */
	      k=idxdval(fx,nxyp,(double)pgx); nxyp++;
	      /* Reallocate memory to increase node array by 1 */
	      if (!(fx=(double *)realloc(fx,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fx\n\
\tarray of size %d",nxyp);
	      if (!(fy=(double *)realloc(fy,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fy\n\
\tarray of size %d",nxyp);
	      if (!(fy2=(double *)realloc(fy2,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fy2\n\
\tarray of size %d",nxyp);
	      /* Insert new node */
	      for (j=nxyp-1; j>k; j--) { fx[j]=fx[j-1]; fy[j]=fy[j-1]; }
	      fx[k]=(double)pgx; fy[k]=(double)pgy;
	    }
	    /* Exit selection loop */
	  }
	  if (!strncmp(pgch,"a",1) ||
	      (par->replay && rp->exit==3 && !strncmp(pgch,"q",1))) {
	    /* Allocate memory for xy-pair array and fill it */
	    if (!(xyp=(twoxyp *)malloc((size_t)(nxyp*sizeof(twoxyp)))))
	      errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\txyp array of size %d",nxyp);
	    for (j=0; j<nxyp; j++) { xyp[j].x1=fx[j]; xyp[j].y1=fy[j]; }
	    /* Remember action by filling in an action report */
	    if (par->replay && rp->exit==3) {
	      (*act)[rp->cact].xyp=xyp; xyp=NULL; (*act)[rp->cact].nxyp=nxyp;
	      (*act)[rp->cact].d[0]=cspec->wl[fsp];
	      (*act)[rp->cact].d[1]=cspec->wl[fep];
	      (*act)[rp->cact].d[2]=cspec->wl[fbsp];
	      (*act)[rp->cact].d[3]=cspec->wl[fbep];
	      (*act)[rp->cact].i[0]=con[i][0]; (*act)[rp->cact].i[1]=con[i][1];
	      nxyp=0;
	    } else {
	      if (!(*nact)) {
		if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
	      }
	      else
		if (!(*act=(action *)
		      realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
	      /* If in paused replay mode then shift later actions to end */
	      if (par->replay && rp->exit==2) {
		for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
		cact=rp->cact++; (*nact)++;
	      } else cact=(*nact)++;
	      (*act)[cact].act=IOACT; (*act)[cact].xyp=xyp; xyp=NULL;
	      (*act)[cact].nxyp=nxyp; nxyp=0; (*act)[cact].d[0]=cspec->wl[fsp];
	      (*act)[cact].d[1]=cspec->wl[fep]; (*act)[cact].d[2]=cspec->wl[fbsp];
	      (*act)[cact].d[3]=cspec->wl[fbep]; (*act)[cact].i[0]=con[i][0];
	      (*act)[cact].i[1]=con[i][1]; (*act)[cact].rcmb=0;
	      (*act)[cact].nordact=1; (*act)[cact].val=1;
	      /* Alter pixel arrays */
	      if (par->thar<=1) {
		for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
		     j<con[i][4]; j++,k++) {
		  scale=spec[con[i][0]].or[con[i][1]].rdco[k]/fco[j];
		  spec[con[i][0]].or[con[i][1]].ordco[k]=scale;
		  spec[con[i][0]].or[con[i][1]].rdco[k]=fco[j];
		  spec[con[i][0]].or[con[i][1]].rdfl[k]*=scale;
		  spec[con[i][0]].or[con[i][1]].rder[k]*=scale;
		  spec[con[i][0]].or[con[i][1]].rdef[k]*=scale;
		  spec[con[i][0]].or[con[i][1]].rdme[k]*=scale;
		}
	      }
	      else {
		for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
		     j<con[i][4]; j++,k++) {
		  spec[con[i][0]].or[con[i][1]].ordco[k]=
		    spec[con[i][0]].or[con[i][1]].rdco[k];
		  spec[con[i][0]].or[con[i][1]].rdco[k]=fco[j];
		}
	      }
	      cp->rccsp=1;
	    }
	  }
	  con_redo=0;
	  /* Clean up temporary arrays */
	  free(fco); free(fx); free(fy); free(fy2);
	  if (par->replay && rp->exit==3) { rp->exit=0; break; }
	}
	if (!plval) {	
	  /** Interpolate new continuum to combined spectrum **/
	  /* Find initial fit region */
	  if (par->replay && rp->exit==3) {
	    /* If in replay mode we must use the values stores in the action */
	    if ((fsp=idxdval(&(cspec->rwl[sp]),cnp,(*act)[rp->cact].d[0]))==-1) {
	      nferrormsg("UVES_plot_cspec(): Cannot find pixel with\n\
\twavelength %lf in combined spectrum",(*act)[rp->cact].d[0]); return 0;
	    } else fsp+=sp;
	    if ((fep=idxdval(&(cspec->rwl[fsp]),cnp+sp-fsp,(*act)[rp->cact].d[1]))
		==-1) fep=ep;
	    else fep+=fsp;
	    if ((fbsp=idxdval(&(cspec->rwl[fsp]),cnp+sp-fsp,(*act)[rp->cact].d[2]))
		==-1) fbsp=fep-2;
	    else fbsp+=fsp;
	    if ((fbep=idxdval(&(cspec->rwl[fbsp]),cnp+sp-fbsp,
			      (*act)[rp->cact].d[3]))==-1) fbep=fep-1;
	    else fbep+=fbsp;
	    fnp=fep-fsp+1;
	    if (fnp<4) {
	      nferrormsg("UVES_plot_cspec(): Too few (=%d) pixels\n\
\tto interpolate across"); return 0;
	    }
	    fbsp=(MAX(fbsp,(fsp+1))); fbep=(MAX(fbep,(fbsp+1)));
	    fep=(MAX(fep,(fbep+1))); fbsnp=fbsp-fsp+1; fbenp=fep-fbep+1;
	    /* Initialize number of x-y pairs */
	    nxyp=(*act)[rp->cact].nxyp;
	  } else {
	    fsp=sp+(int)(FFRAC*cnp); fep=ep-(int)(FFRAC*cnp); fnp=fep-fsp+1;
	    fbsp=fsp+(int)(FBFRAC*(double)fnp); fbsnp=fbsp-fsp+1;
	    fbep=fep-(int)(FBFRAC*(double)fnp); fbenp=fep-fbep+1;
	    nxyp=2;
	  }
          /* Allocate memory for fitting arrays */
          if ((fco=darray(cnp))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fco array of size %d",cnp);
          if ((fx=darray(nxyp))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fx array of size %d",nxyp);
          if ((fy=darray(nxyp))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fy array of size %d",nxyp);
          if ((fy2=darray(nxyp))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fy2 array of size %d",nxyp);
	  if (par->replay && rp->exit==3) {
	    /* Copy over node information and free x-y pair array in
	       original action */
	    for (i=0; i<nxyp; i++) {
	      fx[i]=(*act)[rp->cact].xyp[i].x1; fy[i]=(*act)[rp->cact].xyp[i].y1;
	    }
	    free((*act)[rp->cact].xyp); (*act)[rp->cact].xyp=NULL;
	  } else {
	    /* Initialize y end points of node array to reasonable values */
	    fy[0]=cspec->co[fsp]; fy[1]=cspec->co[fep];
	  }
	  /* Start fit selection loop */
	  /* Fill in fitted continuum array */
	  for (i=0,j=sp; i<cnp; i++,j++) fco[i]=cspec->co[j];
	  /* Start fit selection loop */
	  strcpy(pgch,"\0");
	  while (strncmp(pgch,"q",1) && strncmp(pgch,"a",1)) {
	    /* End points of spline node array are always ends of fitting region */
	    fx[0]=cspec->wl[fsp]; fx[nxyp-1]=cspec->wl[fep];
	    /* Plot raw data again */
	    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd3,cp->vpu3); cpgsci(0);
	    cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
	    cpgsci(15);
	    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
	    cpgsci(6); cpgslw(2.0*plenv.lw);
	    cpgsch(0.9*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
	    cpgsch(0.8*plenv.ch); cpgslw(plenv.lw); cpgsci(3);
	    cpglab(" ",plenv.ylab[3]," "); cpgsci(1);
	    cpgbin(cnp,&(wl[sp]),&(fl[sp]),1);
	    cpgsci(5); cpgsls(4);
	    cpgmove(plenv.xmin[1],0.0); cpgdraw(plenv.xmax[1],0.0); cpgsls(1);
	    /* Reinitialize continuum array */
	    for (i=0,j=sp; i<cnp; i++,j++) fco[i]=cspec->co[j];
	    /* Calculate spline */
	    if (!spline(fx,fy,nxyp,0.0,0.0,fy2,1))
	      errormsg("UVES_plot_cspec(): Error returned from spline()");
	    for (i=fsp,j=fsp-sp; i<=fep; i++,j++)
	      fco[j]=splint(fx,fy,fy2,nxyp,cspec->wl[i]);
	    /* Blend two continua together at ends over the blend region */
	    for (i=0,j=fsp-sp,k=fsp; i<fbsnp; i++,j++,k++) {
	      frac=(double)(i+1)/(double)fbsnp;
	      fco[j]=frac*fco[j]+(1.0-frac)*cspec->co[k];
	    }
	    for (i=0,j=fep-sp,k=fep; i<fbenp; i++,j--,k--) {
	      frac=(double)(i+1)/(double)fbenp;
	      fco[j]=frac*fco[j]+(1.0-frac)*cspec->co[k];
	    }
	    /* Plot the original continuum over the data */
	    for (i=0,j=sp; i<cnp; i++,j++) co[j]=cspec->co[j];
	    cpgsci(3); cpgline(cnp,&(wl[sp]),&(co[sp]));
	    /* Plot new continuum */
	    for (i=fsp-sp,j=fsp; i<fep-sp; i++,j++) co[j]=fco[i];
	    cpgsci(6); cpgline(fnp,&(wl[fsp]),&(co[fsp]));
	    /* Mark spline nodes */
	    cpgsch(0.6*plenv.ch);
	    for (j=0; j<nxyp; j++) cpgpt1((float)fx[j],(float)fy[j],8);
	    /* Plot fit and blend region bounds */
	    cpgsci(6); cpgsls(2);
	    cpgmove(wl[fbsp],plenv.ymin[3]); cpgdraw(wl[fbsp],plenv.ymax[3]);
	    cpgmove(wl[fbep],plenv.ymin[3]); cpgdraw(wl[fbep],plenv.ymax[3]);
	    cpgsls(1);
	    cpgmove(wl[fsp],plenv.ymin[3]); cpgdraw(wl[fsp],plenv.ymax[3]);
	    cpgmove(wl[fep],plenv.ymin[3]); cpgdraw(wl[fep],plenv.ymax[3]);
	    /* Alter fit by selecting new bounds, changing order or
	       selecting/deselecting points */
	    cpgsci(2); cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	    if (!strncmp(pgch,"?",1)) {
	      fprintf(stderr,"\
Options for windowing new plotting region:\n\
 Left mouse  : Select left, right, or smoothing limit to change.\n\
                 Any key to place limit at new position and 'r'\n\
                 registers selected smoothing limit to corresponding\n\
                 fitting limit.\n\
 Right mouse : Remove spline node nearest cursor\n\
 Middle mouse: Move spline node nearest cursor.\n\
                 Any key to place at new position.\n\
 Space       : New spline node at cursor\n\
 a           : Accept fit\n\
 q           : Quit fitting procedure\n\n");
	      cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	    }
	    if (!strncmp(pgch,"A",1)) {
	      clickx=0.01*(plenv.xmax[1]-plenv.xmin[1]);
	      clicknp=idxfval(&(wl[sp]),cnp,wl[sp]+clickx)+2;
	      if (pgx>wl[fsp]-clickx && pgx<wl[fsp]+clickx) {
		cpgsci(15);
		cpgmove(wl[fsp],plenv.ymin[3]); cpgdraw(wl[fsp],plenv.ymax[3]);
		cpgsci(6); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		k=idxfval(&(wl[fsp]),fnp,(float)fx[1])-1; k+=fsp;
		if (!strncmp(pgch,"r",1)) {
		  fsp=(MIN(k,(fbsp-1))); fnp=fep-fsp+1; fbsnp=1;
		}
		else {
		  fsp=((fsp=idxfval(&(wl[sp]),cnp,pgx))==-1) ? fbep-clicknp :
		    fsp+sp; fsp=MIN(fsp,(fbep-clicknp)); fsp=(MIN(k,fsp));
		    fnp=fep-fsp+1;
		    cpgsci(15); cpgmove(wl[fbsp],plenv.ymin[3]);
		    cpgdraw(wl[fbsp],plenv.ymax[3]);
		    fbsp=MAX(fbsp,(MIN((fbep-1),(fsp+fbsnp-1)))); fbsnp=fbsp-fsp+1;
		}
		fy[0]=cspec->co[fsp];
	      }
	      else if (pgx>wl[fep]-clickx && pgx<wl[fep]+clickx) {
		cpgsci(15);
		cpgmove(wl[fep],plenv.ymin[3]); cpgdraw(wl[fep],plenv.ymax[3]);
		cpgsci(6); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		k=idxfval(&(wl[fsp]),fnp,(float)fx[nxyp-2]); k+=fsp;
		if (!strncmp(pgch,"r",1)) {
		  fep=(MAX(k,(fbep+1))); fnp=fep-fsp+1; fbenp=1;
		}
		else {
		  fep=((fep=idxfval(&(wl[sp]),cnp,pgx))==-1) ? ep : fep+sp;
		  fep=MAX(fep,(fbsp+clicknp)); fep=(MAX(k,fep)); fnp=fep-fsp+1;
		  cpgsci(15); cpgmove(wl[fbep],plenv.ymin[3]);
		  cpgdraw(wl[fbep],plenv.ymax[3]);
		  fbep=MIN(fbep,(MAX((fbsp+1),(fep-fbenp+1)))); fbenp=fep-fbep+1;
		}
		fy[nxyp-1]=cspec->co[fep];
	      }
	      else if (pgx>wl[fbsp]-clickx && pgx<wl[fbsp]+clickx) {
		cpgsci(15); cpgmove(wl[fbsp],plenv.ymin[3]);
		cpgdraw(wl[fbsp],plenv.ymax[3]);
		cpgsci(6); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fbsp=fsp+1; fbsnp=1; }
		else {
		  fbsp=((fbsp=idxfval(&(wl[sp]),cnp,pgx))==-1) ? fbep-1: fbsp+sp;
		  fbsp=MAX((fsp+1),(MIN(fbsp,(fbep-1)))); fbsnp=fbsp-fsp+1;
		}
	      }
	      else if (pgx>wl[fbep]-clickx && pgx<wl[fbep]+clickx) {
		cpgsci(15); cpgmove(wl[fbep],plenv.ymin[3]);
		cpgdraw(wl[fbep],plenv.ymax[3]);
		cpgsci(6); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fbep=fep-1; fbenp=1; }
		else {
		  fbep=((fbep=idxfval(&(wl[sp]),cnp,pgx))==-1) ? fep-1 : fbep+sp;
		  fbep=MIN((fep-1),(MAX(fbep,(fbsp+1)))); fbenp=fep-fbep+1;
		}
	      }
	    }
	    else if (!strncmp(pgch,"X",1) && nxyp>2) {
	      /* Remove spline node nearest cursor */
	      k=1; cdist=(pgx-(float)fx[k])/(plenv.xmax[1]-plenv.xmin[1]);
	      cdist=sqrt(cdist*cdist);
	      for (j=2; j<nxyp-1; j++) {
		ftemp=(pgx-(float)fx[j])/(plenv.xmax[1]-plenv.xmin[1]);
		ftemp=sqrt(ftemp*ftemp); if (ftemp<cdist) { cdist=ftemp; k=j; }
	      }
	      nxyp--; for (j=k; j<nxyp; j++) { fx[j]=fx[j+1]; fy[j]=fy[j+1]; }
	      /* Reallocate memory in accordance with removal */
	      if (!(fx=(double *)realloc(fx,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fx\n\
\tarray of size %d",nxyp);
	      if (!(fy=(double *)realloc(fy,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fy\n\
\tarray of size %d",nxyp);
	      if (!(fy2=(double *)realloc(fy2,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fy2\n\
\tarray of size %d",nxyp);
	    }
	    else if (!strncmp(pgch,"D",1)) {
	      /* Move spline node nearest cursor */
	      k=0; cdist=(pgx-(float)fx[k])/(plenv.xmax[1]-plenv.xmin[1]);
	      cdist=sqrt(cdist*cdist);
	      for (j=1; j<nxyp; j++) {
		ftemp=(pgx-(float)fx[j])/(plenv.xmax[1]-plenv.xmin[1]);
		ftemp=sqrt(ftemp*ftemp); if (ftemp<cdist) { cdist=ftemp; k=j; }
	      }
	      /* Erase symbol for selected node */
	      cpgsci(15); cpgsch(0.6*plenv.ch);
	      cpgpt1((float)fx[k],(float)fy[k],8);
	      cpgsci(5); cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	      if (k!=0 && k!=nxyp-1)
		fx[k]=(MAX((MIN(((double)pgx),((1.0-1.e-7)*fx[k+1]))),
			   ((1.0+1.e-7)*fx[k-1])));
	      fy[k]=(double)pgy;
	    }
	    else if (!strncmp(pgch," ",1) && pgx>(float)fx[0] &&
		     pgx<(float)fx[nxyp-1]) {
	      /* Find index at which to insert new node */
	      k=idxdval(fx,nxyp,(double)pgx); nxyp++;
	      /* Reallocate memory to increase node array by 1 */
	      if (!(fx=(double *)realloc(fx,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fx\n\
\tarray of size %d",nxyp);
	      if (!(fy=(double *)realloc(fy,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fy\n\
\tarray of size %d",nxyp);
	      if (!(fy2=(double *)realloc(fy2,(size_t)(nxyp*sizeof(double)))))
		errormsg("UVES_plot_cspec(): Cannot reallocate memory to fy2\n\
\tarray of size %d",nxyp);
	      /* Insert new node */
	      for (j=nxyp-1; j>k; j--) { fx[j]=fx[j-1]; fy[j]=fy[j-1]; }
	      fx[k]=(double)pgx; fy[k]=(double)pgy;
	    }
	    /* Exit selection loop */
	  }
	  if (!strncmp(pgch,"a",1) ||
	      (par->replay && rp->exit==3 && !strncmp(pgch,"q",1))) {
	    /* Allocate memory for xy-pair array and fill it */
	    if (!(xyp=(twoxyp *)malloc((size_t)(nxyp*sizeof(twoxyp)))))
	      errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\txyp array of size %d",nxyp);
	    for (j=0; j<nxyp; j++) { xyp[j].x1=fx[j]; xyp[j].y1=fy[j]; }
	    /* Remember action by filling in an action report */
	    if (par->replay && rp->exit==3) {
	      (*act)[rp->cact].xyp=xyp; xyp=NULL; (*act)[rp->cact].nxyp=nxyp;
	      (*act)[rp->cact].d[0]=cspec->wl[fsp];
	      (*act)[rp->cact].d[1]=cspec->wl[fep];
	      (*act)[rp->cact].d[2]=cspec->wl[fbsp];
	      (*act)[rp->cact].d[3]=cspec->wl[fbep]; nxyp=0;
	    } else {
	      if (!(*nact)) {
		if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
	      }
	      else
		if (!(*act=(action *)
		      realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
	      /* If in paused replay mode then shift later actions to end */
	      if (par->replay && rp->exit==2) {
		for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
		cact=rp->cact++; (*nact)++;
	      } else cact=(*nact)++;
	      (*act)[cact].act=ICACT; (*act)[cact].xyp=xyp; xyp=NULL;
	      (*act)[cact].nxyp=nxyp; nxyp=0; (*act)[cact].d[0]=cspec->wl[fsp];
	      (*act)[cact].d[1]=cspec->wl[fep]; (*act)[cact].d[2]=cspec->wl[fbsp];
	      (*act)[cact].d[3]=cspec->wl[fbep]; (*act)[cact].rcmb=0;
	      (*act)[cact].nordact=1; (*act)[cact].val=1;
	      /* Alter relevant pixel values */
	      if (par->thar<=1) {
		for (i=0,j=sp; i<cnp; i++,j++) {
		  cspec->oco[j]=cspec->co[j]; scale=cspec->co[j]/fco[i];
		  fscale=(float)scale; cspec->co[j]=fco[i]; cspec->no[j]*=scale;
		  cspec->ne[j]*=scale; cspec->nf[j]*=scale; co[j]=fco[i];
		  no[j]*=fscale; ne[j]*=fscale;
		  for (k=0; k<ncon; k++) {
		    if (con[k][2]<=j && con[k][3]>=j) {
		      idx=j-spec[con[k][0]].or[con[k][1]].csidx;
		      spec[con[k][0]].or[con[k][1]].rdco[idx]=cspec->co[j];
		    }
		  }
		}
	      } else {
		for (i=0,j=sp; i<cnp; i++,j++) {
		  cspec->oco[j]=cspec->co[j]; scale=cspec->co[j]-fco[i];
		  fscale=(float)scale; cspec->co[j]=fco[i]; cspec->no[j]+=scale;
		  co[j]=fco[i]; no[j]+=fscale;
		  for (k=0; k<ncon; k++) {
		    if (con[k][2]<=j && con[k][3]>=j) {
		      idx=j-spec[con[k][0]].or[con[k][1]].csidx;
		      spec[con[k][0]].or[con[k][1]].rdco[idx]=cspec->co[j];
		    }
		  }
		}
	      }
	    }
	  }
	  con_redo=0;
	  /* Clean up temporary arrays */
	  free(fco); free(fx); free(fy); free(fy2);
	  if (par->replay && rp->exit==3) { rp->exit=0; break; }
	  if (!strncmp(pgch,"q",1)) for (i=fsp; i<=fep; i++) co[i]=cspec->co[i];
	}
      }
    } else if (!strncmp(pgch,"f",1) ||
	       (!strncmp(pgch,"m",1) && plval && npcon==1 && cact>-1 &&
		(*act)[cact].act==FOACT)) {
      /*** Fit a new continuum, either to the final combined spectrum or to
	  a single order (if only one plotted) ***/
      /* Decide whether fitting is feasible given current plot state */
      if (!plval || npcon==1) {
	/* Refresh number of pairs of xy pairs of positions */
	if (plval && npcon==1) {
	/** Fit new continuum for a single order **/
	  /* Find which spectrum and order we're dealing with */
	  i=0; while (!con[i][6]) i++;
          /* Allocate memory for fitting arrays */
          if ((fco=darray(con[i][4]))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fco array of size %d",con[i][4]);
          if ((fst=iarray(con[i][4]))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fst array of size %d",con[i][4]);
	  /* Find initial fit region */
	  m=0; /* Counter which flags mirroring of previous fit on current order */
	  if (par->replay && rp->exit==3) {
	    /* If in replay mode we must use the values stored in the action */
	    if ((fsp=idxdval(&(cspec->rwl[sp]),cnp,(*act)[rp->cact].d[0]))==-1) {
	      nferrormsg("UVES_plot_cspec(): Cannot find pixel with\n\
\twavelength %lf in combined spectrum",(*act)[rp->cact].d[0]); return 0;
	    } else fsp+=sp;
	    if ((fep=idxdval(&(cspec->rwl[fsp]),cnp+sp-fsp,(*act)[rp->cact].d[1]))
		==-1) fep=ep;
	    else fep+=fsp;
	    if ((fbsp=idxdval(&(cspec->rwl[fsp]),cnp+sp-fsp,(*act)[rp->cact].d[2]))
		==-1) fbsp=fep-2;
	    else fbsp+=fsp;
	    if ((fbep=idxdval(&(cspec->rwl[fbsp]),cnp+sp-fbsp,
			      (*act)[rp->cact].d[3]))==-1) fbep=fep-1;
	    else fbep+=fbsp;
	    fnp=fep-fsp+1;
	    if (fnp<4) {
	      nferrormsg("UVES_plot_cspec(): Too few (=%d) pixels\n\
\tto be fit"); return 0;
	    }
	    fbsp=(MAX(fbsp,(fsp+1))); fbep=(MAX(fbep,(fbsp+1)));
	    fep=(MAX(fep,(fbep+1))); fbsnp=fbsp-fsp+1; fbenp=fep-fbep+1;
	    rejsigl=(*act)[rp->cact].d[4]; rejsigu=(*act)[rp->cact].d[5];
	    fit_typ=(*act)[rp->cact].i[2]; fit_ord=(*act)[rp->cact].i[3];
	    /* Direct memory for x-y pair arrays */
	    nxyp=(*act)[rp->cact].nxyp; xyp=(*act)[rp->cact].xyp;
	    /* Initialize status array with rejected or included pixels */
	    for (j=0,k=con[i][2]; j<con[i][4]; j++,k++) {
	      for (l=0; l<nxyp; l++) {
		if (wl[k]>xyp[l].x1 && wl[k]<xyp[l].x2 && fl[k]>xyp[l].y1 &&
		    fl[k]<xyp[l].y2) fst[j]=(xyp[l].i) ? 1 : PCLIP;
	      }
	    }
	  } else if (!strncmp(pgch,"m",1) && cact>-1 && (*act)[cact].act==FOACT) {
	    /* If mirroring previous action, i.e. a fit to an order,
	       check that that action's parameters make sense for the
	       current order */
	    if (cspec->wl[sp]<(*act)[cact].d[0] && cspec->wl[ep]>(*act)[cact].d[1]) {
	      /* Determine whether previous fit started at order edge -
		 mirror that situation for the current order. If not,
		 start current fit at same wavelength as previous one. */
	      wl1=cspec->wl[spec[(*act)[cact].i[0]].or[(*act)[cact].i[1]].csidx];
	      wl2=(par->linear) ? wl1+0.5*par->disp : wl1*(1.0+0.5*par->disp/C_C_K);
	      wl1=(par->linear) ? wl1-0.5*par->disp : wl1*(1.0-0.5*par->disp/C_C_K);
	      if ((*act)[cact].d[0]>wl1 && (*act)[cact].d[0]<wl2) {
		if (con[i][2]==spec[con[i][0]].or[con[i][1]].csidx) fsp=con[i][2];
		else strcpy(pgch,"q\0");
	      } else {
		if ((fsp=idxdval(&(cspec->wl[con[i][2]]),con[i][4],(*act)[cact].d[0]))
		    ==-1) strcpy(pgch,"q\0");
		else fsp+=con[i][2];
	      }
	      /* Now determine most suitable place to end left-hand
		 blending region */
	      wl2=(par->linear) ? wl2+par->disp : wl2*(1.0+par->disp/C_C_K);
	      if ((*act)[cact].d[2]>wl1 && (*act)[cact].d[2]<wl2) fbsp=fsp+1;
	      else if (fsp>=con[i][2] && fsp<con[i][3]) {
		if ((fbsp=idxdval(&(cspec->wl[fsp]),con[i][3]-fsp+1,(*act)[cact].d[2]))
		    ==-1) strcpy(pgch,"q\0");
		else { fbsp+=fsp; fbsp=(MAX(fsp+1,fbsp)); }
	      } else strcpy(pgch,"q\0");
	      /* Determine whether previous fit ended at order edge -
		 mirror that situation for the current order. If not,
		 end current fit at same wavelength as previous one. */
	      wl1=cspec->wl[spec[(*act)[cact].i[0]].or[(*act)[cact].i[1]].ceidx];
	      wl2=(par->linear) ? wl1+0.5*par->disp : wl1*(1.0+0.5*par->disp/C_C_K);
	      wl1=(par->linear) ? wl1-0.5*par->disp : wl1*(1.0-0.5*par->disp/C_C_K);
	      if ((*act)[cact].d[1]>wl1 && (*act)[cact].d[1]<wl2) {
		if (con[i][3]==spec[con[i][0]].or[con[i][1]].ceidx) fep=con[i][3];
		else strcpy(pgch,"q\0");
	      } else {
		if ((fep=idxdval(&(cspec->wl[con[i][2]]),con[i][4],(*act)[cact].d[1]))
		    ==-1) fep=con[i][3];
		else fep+=con[i][2];
	      }
	      /* Now determine most suitable place to begin right-hand
		 blending region */
	      wl1=(par->linear) ? wl1-par->disp : wl1*(1.0-par->disp/C_C_K);
	      if ((*act)[cact].d[3]>wl1 && (*act)[cact].d[3]<wl2) fbep=fep-1;
	      else if (fep>con[i][2] && fep<=con[i][3]) {
		if ((fbep=idxdval(&(cspec->wl[fsp]),con[i][3]-fsp+1,(*act)[cact].d[3]))
		    ==-1) fbep=fep-1;
		else { fbep+=fsp; fbep=(MIN(fep-1,fbep)); }
	      } else strcpy(pgch,"q\0");
	      fnp=fep-fsp+1;
	      if (fbsp>=fep || fbsp>fbep || fbep<=fsp || fbep<fbsp) strcpy(pgch,"q\0");
	      fbsnp=fbsp-fsp+1; fbenp=fep-fbep+1;
	      fit_typ=(*act)[cact].i[2]; fit_ord=(*act)[cact].i[3];
	      rejsigu=(*act)[cact].d[5]; rejsigl=(*act)[cact].d[4];
	      /* Indicates to fit selection loop to consider
		 (de)selecting points from fit */
	      m=(*act)[cact].nxyp;
	    } else strcpy(pgch,"q\0");
	  } else {
	    fsp=(con[i][2]==spec[con[i][0]].or[con[i][1]].csidx) ? con[i][2] :
	      con[i][2]+(int)(FFRAC*con[i][4]);
	    fep=(con[i][3]==spec[con[i][0]].or[con[i][1]].ceidx) ? con[i][3] :
	      con[i][3]-(int)(FFRAC*con[i][4]);
	    fnp=fep-fsp+1;
	    fbsp=fsp+(int)(FBFRAC*(double)fnp); fbsnp=fbsp-fsp+1;
	    fbep=fep-(int)(FBFRAC*(double)fnp); fbenp=fep-fbep+1;
	    fit_typ=par->contftyp; fit_ord=(MIN(par->contord,fnp));
	    rejsigu=par->contsigu; rejsigl=par->contsigl;
	    nxyp=0;
	  }
	  /* Start fit selection loop */
	  if (strncmp(pgch,"q",1)) strcpy(pgch,"\0");
	  while (strncmp(pgch,"q",1) && strncmp(pgch,"a",1)) {
	    fosp=fsp-spec[con[i][0]].or[con[i][1]].csidx; foep=fosp+fnp-1;
	    /* Plot raw data again */
	    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd3,cp->vpu3); cpgsci(0);
	    cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
	    cpgsci(15);
	    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
	    cpgsci(4); cpgslw(2.0*plenv.lw);
	    cpgsch(0.9*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
	    cpgsch(0.8*plenv.ch); cpgslw(plenv.lw); cpgsci(3);
	    cpglab(" ",plenv.ylab[3]," "); cpgsci(con[i][5]);
	    for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
		 j<con[i][4]; j++,k++)
	      or[j]=spec[con[i][0]].or[con[i][1]].rdfl[k];
	    cpgbin(con[i][4],&(wl[con[i][2]]),or,1);
	    /* Fill status array and plot status, including symbols for
	       deselected points */
	    for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx,l=con[i][2];
		 j<con[i][4]; j++,k++,l++) {
	      fco[j]=spec[con[i][0]].or[con[i][1]].rdco[k];
	      if (!fst[j])
		fst[j]=(spec[con[i][0]].or[con[i][1]].rdst[k]==SCLIP ||
			spec[con[i][0]].or[con[i][1]].rdst[k]==LCLIP)
		  ? 1 : spec[con[i][0]].or[con[i][1]].rdst[k];
	      switch (fst[j]) {
	      case RCLIP: cpgpt(1,&(wl[l]),&(or[j]),12); break;
	      case OCLIP: cpgpt(1,&(wl[l]),&(or[j]),3); break;
	      case ACLIP: cpgpt(1,&(wl[l]),&(or[j]),8); break;
	      case NCLIP: cpgpt(1,&(wl[l]),&(or[j]),5); break;
	      case CCLIP: cpgpt(1,&(wl[l]),&(or[j]),7); break;
	      case PCLIP: cpgpt(1,&(wl[l]),&(or[j]),14); break;
	      }
	    }
	    cpgsci(5); cpgsls(4);
	    cpgmove(plenv.xmin[1],0.0); cpgdraw(plenv.xmax[1],0.0); cpgsls(1);
	    /* Do a fit for the selected data */
	    fit_rs=UVES_confit(&(spec[con[i][0]].or[con[i][1]].rdfl[fosp]),
			       &(spec[con[i][0]].or[con[i][1]].rder[fosp]),
			       &(fst[fsp-con[i][2]]),fnp,fit_typ,fit_ord,rejsigl,
			       rejsigu,0.0,1,&(fco[fsp-con[i][2]]));
	    if (fit_rs==-1 || fit_rs==0) { strcpy(pgch,"q\0"); break; }
	    /* Blend two continua together at ends over the blend region */
	    for (j=0,k=fsp-con[i][2],l=fosp; j<fbsnp; j++,k++,l++) {
	      frac=(double)(j+1)/(double)fbsnp;
	      fco[k]=frac*fco[k]+(1.0-frac)*spec[con[i][0]].or[con[i][1]].rdco[l];
	    }
	    for (j=0,k=fep-con[i][2],l=foep; j<fbenp; j++,k--,l--) {
	      frac=(double)(j+1)/(double)fbenp;
	      fco[k]=frac*fco[k]+(1.0-frac)*spec[con[i][0]].or[con[i][1]].rdco[l];
	    }
	    /* Plot the original continuum over the data */
	    cpgsci((PMOD(con[i][5],1,14,2))); cpgline(cnp,&(wl[sp]),&(co[sp]));
	    /* Plot new continuum */
	    for (j=0,k=fsp-con[i][2]; j<fnp; j++,k++) or[j]=fco[k];
	    cpgsci((PMOD(con[i][5],2,14,2))); cpgline(fnp,&(wl[fsp]),or);
	    /* Plot fit and blend region bounds */
	    cpgsci(3); cpgsls(2);
	    cpgmove(wl[fbsp],plenv.ymin[3]); cpgdraw(wl[fbsp],plenv.ymax[3]);
	    cpgmove(wl[fbep],plenv.ymin[3]); cpgdraw(wl[fbep],plenv.ymax[3]);
	    cpgsls(1);
	    cpgmove(wl[fsp],plenv.ymin[3]); cpgdraw(wl[fsp],plenv.ymax[3]);
	    cpgmove(wl[fep],plenv.ymin[3]); cpgdraw(wl[fep],plenv.ymax[3]);
	    /* Alter fit by selecting new bounds, changing fit type, order or
	       selecting/deselecting points */
	    cpgsci((PMOD(con[i][5],3,14,2)));
	    if (m) {
	      pgx=(*act)[cact].xyp[m-1].x1; pgy=(*act)[cact].xyp[m-1].y1;
	      if ((*act)[cact].xyp[m-1].i) strcpy(pgch,"D\0");
	      else strcpy(pgch,"X\0");
	    } else cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	    if (!strncmp(pgch,"?",1)) {
	      fprintf(stderr,"\
Options for windowing new plotting region:\n\
 Left mouse  : Select left, right, or smoothing limit to change.\n\
                 Any key to place limit at new position and 'r'\n\
                 registers selected smoothing limit to corresponding\n\
                 fitting limit.\n\
 Right mouse : Clip pixels from subsequent continuum fit.\n\
                 Right mouse again selects second corner of clip window,\n\
                 any other entry to abort selection.\n\
 Middle mouse: Unclip pixels from subsequent continuum fit.\n\
                 Middle mouse again selects second corner of clip window,\n\
                 any other entry to abort selection.\n\
 +           : Increase fit order by 1.\n\
 -           : Decrease fit order by 1.\n\
 a           : Accept fit\n\
 c           : Cancel any and all clip and unclip regions\n\
 l, u        : Alter lower/upper sigma-clipping rejection threshold.\n\
                 Points more than %5.2lf-sigma below & %5.2lf-sigma above curve\n\
                 are rejected at each iteration.\n\
 t           : Alter fit type. 1=polynomial, 2=Chebyshev, 3=Legendre.\n\
                 Current fit type is %d.\n\
 o           : Alter fit order. Current fit order is %d.\n\
 q           : Quit fitting procedure\n\n",rejsigl,rejsigu,fit_typ,fit_ord);
	      cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	    }
	    if (!strncmp(pgch,"A",1)) {
	      clickx=0.01*(plenv.xmax[1]-plenv.xmin[1]);
	      clicknp=idxfval(&(wl[con[i][2]]),con[i][4],wl[con[i][2]]+clickx);
	      clicknp+=2;
	      if (pgx>wl[fsp]-clickx && pgx<wl[fsp]+clickx) {
		cpgsci(15);
		cpgmove(wl[fsp],plenv.ymin[3]); cpgdraw(wl[fsp],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fsp=fbsp-1; fbsnp=1; }
		else {
		  fsp=((fsp=idxfval(&(wl[con[i][2]]),con[i][4],pgx))==-1) ?
		    fbep-clicknp : fsp+con[i][2]; fsp=MIN(fsp,(fbep-clicknp));
		  fnp=fep-fsp+1;
		  cpgsci(15); cpgmove(wl[fbsp],plenv.ymin[3]);
		  cpgdraw(wl[fbsp],plenv.ymax[3]);
		  fbsp=MAX(fbsp,(MIN((fbep-1),(fsp+fbsnp-1)))); fbsnp=fbsp-fsp+1;
		}
	      }
	      else if (pgx>wl[fep]-clickx && pgx<wl[fep]+clickx) {
		cpgsci(15);
		cpgmove(wl[fep],plenv.ymin[3]); cpgdraw(wl[fep],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fep=fbep+1; fbenp=1; }
		else {
		  fep=((fep=idxfval(&(wl[con[i][2]]),con[i][4],pgx))==-1) ?
		    con[i][3] : fep+con[i][2]; fep=MAX(fep,(fbsp+clicknp));
		  fnp=fep-fsp+1;
		  cpgsci(15); cpgmove(wl[fbep],plenv.ymin[3]);
		  cpgdraw(wl[fbep],plenv.ymax[3]);
		  fbep=MIN(fbep,(MAX((fbsp+1),(fep-fbenp+1)))); fbenp=fep-fbep+1;
		}
	      }
	      else if (pgx>wl[fbsp]-clickx && pgx<wl[fbsp]+clickx) {
		cpgsci(15); cpgmove(wl[fbsp],plenv.ymin[3]);
		cpgdraw(wl[fbsp],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fbsp=fsp+1; fbsnp=1; }
		else {
		  fbsp=((fbsp=idxfval(&(wl[con[i][2]]),con[i][4],pgx))==-1) ?
		    fbep-1 : fbsp+con[i][2];
		  fbsp=MAX((fsp+1),(MIN(fbsp,(fbep-1)))); fbsnp=fbsp-fsp+1;
		}
	      }
	      else if (pgx>wl[fbep]-clickx && pgx<wl[fbep]+clickx) {
		cpgsci(15); cpgmove(wl[fbep],plenv.ymin[3]);
		cpgdraw(wl[fbep],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fbep=fep-1; fbenp=1; }
		else {
		  fbep=((fbep=idxfval(&(wl[con[i][2]]),con[i][4],pgx))==-1) ?
		    fep-1 : fbep+con[i][2];
		  fbep=MIN((fep-1),(MAX(fbep,(fbsp+1)))); fbenp=fep-fbep+1;
		}
	      }
	    }
	    else if (!strncmp(pgch,"X",1)) {
	      /* Reject pixels from next fit */
	      cpgsci(2); pgxo=pgx; pgyo=pgy;
	      if (m) {
		pgx=(*act)[cact].xyp[m-1].x2; pgy=(*act)[cact].xyp[m-1].y2;
		strcpy(pgch,"X\0"); m--;
	      } else cpgband(2,0,pgx,pgy,&pgx,&pgy,pgch);
	      if (!strncmp(pgch,"X",1)) {
		if (pgxo>pgx) { FSWAP(pgxo,pgx); }
		if (pgyo>pgy) { FSWAP(pgyo,pgy); }
		for (j=0,selpix=0,k=con[i][2],
		       l=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
		     j<con[i][4]; j++,k++,l++) {
		  if (fst[j]==1 && (wl[k]>pgxo && wl[k]<pgx) &&
		      (spec[con[i][0]].or[con[i][1]].rdfl[l]>(double)pgyo &&
		       spec[con[i][0]].or[con[i][1]].rdfl[l]<(double)pgy)) {
		    fst[j]=PCLIP; selpix++;
		  }
		}
		/* Allocate memory for a new pair of xy pairs */
		if (!nxyp && selpix) {
		  if (!(xyp=(twoxyp *)malloc((size_t)(sizeof(twoxyp)))))
		    errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\txyp array of size 1");
		}
		else if (selpix) {
		  if (!(xyp=(twoxyp *)
			realloc(xyp,(size_t)((nxyp+1)*sizeof(twoxyp)))))
		    errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\txyp array to size %d",nxyp+1);
		}
		if (selpix) {
		  xyp[nxyp].x1=pgxo; xyp[nxyp].x2=pgx;
		  xyp[nxyp].y1=pgyo; xyp[nxyp].y2=pgy;
		  xyp[nxyp++].i=0;
		}
	      }
	    }
	    else if (!strncmp(pgch,"D",1)) {
	      /* Include pixels in next fit */
	      cpgsci(2); pgxo=pgx; pgyo=pgy;
	      if (m) {
		pgx=(*act)[cact].xyp[m-1].x2; pgy=(*act)[cact].xyp[m-1].y2;
		strcpy(pgch,"D\0"); m--;
	      } else cpgband(2,0,pgx,pgy,&pgx,&pgy,pgch);
	      if (!strncmp(pgch,"D",1)) {
		if (pgxo>pgx) { FSWAP(pgxo,pgx); }
		if (pgyo>pgy) { FSWAP(pgyo,pgy); }
		for (j=0,selpix=0,k=con[i][2],
		       l=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
		     j<con[i][4]; j++,k++,l++) {
		  if (fst[j]==PCLIP && (wl[k]>pgxo && wl[k]<pgx) &&
		      (spec[con[i][0]].or[con[i][1]].rdfl[l]>(double)pgyo &&
		       spec[con[i][0]].or[con[i][1]].rdfl[l]<(double)pgy)) {
		    fst[j]=1; selpix++;
		  }
		}
		/* Allocate memory for a new pair of xy pairs */
		if (!nxyp && selpix) {
		  if (!(xyp=(twoxyp *)malloc((size_t)(sizeof(twoxyp)))))
		    errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\txyp array of size 1");
		}
		else if (selpix) {
		  if (!(xyp=(twoxyp *)
			realloc(xyp,(size_t)((nxyp+1)*sizeof(twoxyp)))))
		    errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\txyp array to size %d",nxyp+1);
		}
		if (selpix) {
		  xyp[nxyp].x1=pgxo; xyp[nxyp].x2=pgx;
		  xyp[nxyp].y1=pgyo; xyp[nxyp].y2=pgy;
		  xyp[nxyp++].i=1;
		}
	      }
	    }
	    else if (!strncmp(pgch,"c",1)) {
	      for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
		   j<con[i][4]; j++,k++)
		fst[j]=spec[con[i][0]].or[con[i][1]].rdst[k];
	      if (nxyp) { free(xyp); nxyp=0; }
	    }
	    else if (!strncmp(pgch,"t",1)) {
	      firstfit=1; while (firstfit || fit_typ<FITPOL || fit_typ>FITLEG) {
		get_input("Continuum fit type? [1=Poly; 2=Cheb; 3=Lege]","%d",
			  &fit_typ);
		if (fit_typ<FITPOL || fit_typ>FITLEG)
		  warnmsg("UVES_plot_cspec(): Fit type must be >=%d and <=%d",
			  FITPOL,FITLEG);
		firstfit=0;
	      }
	    }
	    else if (!strncmp(pgch,"o",1)) {
	      firstfit=1; while (firstfit || (fit_ord<1 || fnp<=fit_ord)) {
		get_input("Order of continuum fit?","%d",&fit_ord);
		if (fit_ord<=0)
		  warnmsg("UVES_plot_cspec(): Fit order must be > 0");
		if (fnp<=fit_ord)
		  warnmsg("UVES_plot_cspec(): Not enough valid data points\n\
\t(=%d) in selected area to fit an %d-order polynomial",fnp,fit_ord);
		firstfit=0;
	      }
	    }
	    else if (!strncmp(pgch,"+",1)) {
	      fit_ord++;
	      if (fnp<=fit_ord) {
		warnmsg("UVES_plot_cspec(): Not enough valid data points\n\
\t(=%d) in selected area to fit an %d-order polynomial",fnp,fit_ord);
		fit_ord--;
	      }
	    }
	    else if (!strncmp(pgch,"-",1)) {
	      fit_ord--;
	      if (fit_ord<1) {
		warnmsg("UVES_plot_cspec(): Fit order must be > 0");
		fit_ord++;
	      }
	    }
	    else if (!strncmp(pgch,"u",1)) {
	      firstfit=1; while (firstfit || rejsigu<DRNDTOL) {
		get_input("Upper rejection threshold in sigma?","%lf",&rejsigu);
		if (rejsigu<DRNDTOL)
		  warnmsg("UVES_plot_cspec(): Rej.-threshold must be > %5.2lg",
			  DRNDTOL);
		firstfit=0;
	      }
	    }
	    else if (!strncmp(pgch,"l",1)) {
	      firstfit=1; while (firstfit || rejsigl<DRNDTOL) {
		get_input("Lower rejection threshold in sigma?","%lf",&rejsigl);
		if (rejsigl<DRNDTOL)
		  warnmsg("UVES_plot_cspec(): Rej.-threshold must be > %5.2lg",
			  DRNDTOL);
		firstfit=0;
	      }
	    }
	    /* Exit selection loop */
	  }
	  if (!strncmp(pgch,"a",1) ||
	      (par->replay && rp->exit==3 && !strncmp(pgch,"q",1))) {
	    /* Remember action by filling in an action report */
	    if (par->replay && rp->exit==3) {
	      (*act)[rp->cact].xyp=xyp; xyp=NULL;
	      (*act)[rp->cact].nxyp=nxyp; nxyp=0;
	      (*act)[rp->cact].d[0]=cspec->wl[fsp];
	      (*act)[rp->cact].d[1]=cspec->wl[fep];
	      (*act)[rp->cact].d[2]=cspec->wl[fbsp];
	      (*act)[rp->cact].d[3]=cspec->wl[fbep];
	      (*act)[rp->cact].d[4]=rejsigl; (*act)[rp->cact].d[5]=rejsigu;
	      (*act)[rp->cact].i[0]=con[i][0]; (*act)[rp->cact].i[1]=con[i][1];
	      (*act)[rp->cact].i[2]=fit_typ; (*act)[rp->cact].i[3]=fit_ord;
	    } else {
	      if (!(*nact)) {
		if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
	      }
	      else
		if (!(*act=(action *)
		      realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
	      /* If in paused replay mode then shift later actions to end */
	      if (par->replay && rp->exit==2) {
		for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
		cact=rp->cact++; (*nact)++;
	      } else cact=(*nact)++;
	      (*act)[cact].act=FOACT; (*act)[cact].xyp=xyp; xyp=NULL;
	      (*act)[cact].nxyp=nxyp; nxyp=0; (*act)[cact].d[0]=cspec->wl[fsp];
	      (*act)[cact].d[1]=cspec->wl[fep]; (*act)[cact].d[2]=cspec->wl[fbsp];
	      (*act)[cact].d[3]=cspec->wl[fbep]; (*act)[cact].d[4]=rejsigl;
	      (*act)[cact].d[5]=rejsigu; (*act)[cact].i[0]=con[i][0];
	      (*act)[cact].i[1]=con[i][1]; (*act)[cact].i[2]=fit_typ;
	      (*act)[cact].i[3]=fit_ord; (*act)[cact].rcmb=0;
	      (*act)[cact].nordact=1; (*act)[cact].val=1;
	      /* Alter pixel arrays */
	      if (par->thar<=1) {
		for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
		     j<con[i][4]; j++,k++) {
		  scale=spec[con[i][0]].or[con[i][1]].rdco[k]/fco[j];
		  spec[con[i][0]].or[con[i][1]].ordco[k]=scale;
		  spec[con[i][0]].or[con[i][1]].rdco[k]=fco[j];
		  spec[con[i][0]].or[con[i][1]].rdfl[k]*=scale;
		  spec[con[i][0]].or[con[i][1]].rder[k]*=scale;
		  spec[con[i][0]].or[con[i][1]].rdef[k]*=scale;
		  spec[con[i][0]].or[con[i][1]].rdme[k]*=scale;
		}
	      }
	      else {
		for (j=0,k=con[i][2]-spec[con[i][0]].or[con[i][1]].csidx;
		     j<con[i][4]; j++,k++) {
		  spec[con[i][0]].or[con[i][1]].ordco[k]=
		    spec[con[i][0]].or[con[i][1]].rdco[k];
		  spec[con[i][0]].or[con[i][1]].rdco[k]=fco[j];
		}
	      }
	      cp->rccsp=1;
	    }
	  }
	  con_redo=0;
	  /* Clean up temporary arrays */
	  free(fco); free(fst);
	  if (par->replay && rp->exit==3) { rp->exit=0; break; }
	  if (!strncmp(pgch,"q",1)) { if (nxyp) { free(xyp); nxyp=0; } }
	}
	if (!plval) {	
	/** Fit new continuum to combined spectrum **/
          /* Allocate memory for fitting arrays */
          if ((fco=darray(cnp))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fco array of size %d",cnp);
          if ((fst=iarray(cnp))==NULL)
            errormsg("UVES_plot_cspec(): Could not allocate memory\n\
\tfor fst array of size %d",cnp);
	  /* Find initial fit region */
	  if (par->replay && rp->exit==3) {
	    /* If in replay mode we must use the values stores in the action */
	    if ((fsp=idxdval(&(cspec->rwl[sp]),cnp,(*act)[rp->cact].d[0]))==-1) {
	      nferrormsg("UVES_plot_cspec(): Cannot find pixel with\n\
\twavelength %lf in combined spectrum",(*act)[rp->cact].d[0]); return 0;
	    } else fsp+=sp;
	    if ((fep=idxdval(&(cspec->rwl[fsp]),cnp+sp-fsp,(*act)[rp->cact].d[1]))
		==-1) fep=ep;
	    else fep+=fsp;
	    if ((fbsp=idxdval(&(cspec->rwl[fsp]),cnp+sp-fsp,(*act)[rp->cact].d[2]))
		==-1) fbsp=fep-2;
	    else fbsp+=fsp;
	    if ((fbep=idxdval(&(cspec->rwl[fbsp]),cnp+sp-fbsp,
			      (*act)[rp->cact].d[3]))==-1) fbep=fep-1;
	    else fbep+=fbsp;
	    fnp=fep-fsp+1;
	    if (fnp<4) {
	      nferrormsg("UVES_plot_cspec(): Too few (=%d) pixels\n\
\tto be fit"); return 0;
	    }
	    fbsp=(MAX(fbsp,(fsp+1))); fbep=(MAX(fbep,(fbsp+1)));
	    fep=(MAX(fep,(fbep+1))); fbsnp=fbsp-fsp+1; fbenp=fep-fbep+1;
	    rejsigl=(*act)[rp->cact].d[4]; rejsigu=(*act)[rp->cact].d[5];
	    fit_typ=(*act)[rp->cact].i[0]; fit_ord=(*act)[rp->cact].i[1];
	    /* Direct memory for x-y pair arrays */
	    nxyp=(*act)[rp->cact].nxyp; xyp=(*act)[rp->cact].xyp;
	    /* Initialize status array with rejected or included pixels */
	    for (i=0,j=sp; i<cnp; i++,j++) {
	      for (k=0; k<nxyp; k++) {
		if (wl[j]>xyp[k].x1 && wl[j]<xyp[k].x2 && fl[j]>xyp[k].y1 &&
		    fl[j]<xyp[k].y2) fst[i]=(xyp[k].i) ? 1 : PCLIP;
	      }
	    }
	  } else {
	    fsp=sp+(int)(FFRAC*cnp); fep=ep-(int)(FFRAC*cnp); fnp=fep-fsp+1;
	    fbsp=fsp+(int)(FBFRAC*(double)fnp); fbsnp=fbsp-fsp+1;
	    fbep=fep-(int)(FBFRAC*(double)fnp); fbenp=fep-fbep+1;
	    fit_typ=par->contftyp; fit_ord=(MIN(par->contord,fnp));
	    rejsigu=par->contsigu; rejsigl=par->contsigl; nxyp=0;
	  }
	  /* Start fit selection loop */
	  strcpy(pgch,"\0");
	  while (strncmp(pgch,"q",1) && strncmp(pgch,"a",1)) {
	    /* Plot raw data again */
	    cpgsvp(plenv.vpl,plenv.vpr,cp->vpd3,cp->vpu3); cpgsci(0);
	    cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
	    cpgsci(15);
	    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[3],plenv.ymax[3]);
	    cpgsci(4); cpgslw(2.0*plenv.lw);
	    cpgsch(0.9*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
	    cpgsch(0.8*plenv.ch); cpgslw(plenv.lw); cpgsci(3);
	    cpglab(" ",plenv.ylab[3]," "); cpgsci(1);
	    cpgbin(cnp,&(wl[sp]),&(fl[sp]),1);
	    /* Fill status array and plot status, including symbols for
	       deselected points */
	    for (i=0,j=sp; i<cnp; i++,j++) {
	      fco[i]=cspec->co[j];
	      if (!fst[i]) fst[i]=cspec->st[j];
	      switch (fst[i]) {
	      case NCLIP: cpgpt(1,&(wl[j]),&(fl[j]),5); break;
	      case CCLIP: cpgpt(1,&(wl[j]),&(fl[j]),7); break;
	      case PCLIP: cpgpt(1,&(wl[j]),&(fl[j]),14); break;
	      }
	    }
	    cpgsci(5); cpgsls(4);
	    cpgmove(plenv.xmin[1],0.0); cpgdraw(plenv.xmax[1],0.0); cpgsls(1);
	    /* Do a fit for the selected data */
	    fit_rs=UVES_confit(&(cspec->fl[fsp]),&(cspec->er[fsp]),&(fst[fsp-sp]),
			       fnp,fit_typ,fit_ord,rejsigl,rejsigu,0.0,1,
			       &(fco[fsp-sp]));
	    if (fit_rs==-1 || fit_rs==0) { strcpy(pgch,"q\0"); break; }
	    /* Blend two continua together at ends over the blend region */
	    for (i=0,j=fsp-sp,k=fsp; i<fbsnp; i++,j++,k++) {
	      frac=(double)(i+1)/(double)fbsnp;
	      fco[j]=frac*fco[j]+(1.0-frac)*cspec->co[k];
	    }
	    for (i=0,j=fep-sp,k=fep; i<fbenp; i++,j--,k--) {
	      frac=(double)(i+1)/(double)fbenp;
	      fco[j]=frac*fco[j]+(1.0-frac)*cspec->co[k];
	    }
	    /* Fill in rest of continuum array */
	    for (i=0,j=sp; i<fsp-sp; i++,j++) fco[i]=cspec->co[j];
	    for (i=fnp-1,j=ep-1; i>fep; i--,j--) fco[i]=cspec->co[j];
	    /* Plot the original continuum over the data */
	    for (i=0,j=sp; i<cnp; i++,j++) co[j]=cspec->co[j];
	    cpgsci(3); cpgline(cnp,&(wl[sp]),&(co[sp]));
	    /* Plot new continuum */
	    for (i=fsp-sp,j=fsp; i<fep-sp; i++,j++) co[j]=fco[i];
	    cpgsci(4); cpgline(fnp,&(wl[fsp]),&(co[fsp]));
	    /* Plot fit and blend region bounds */
	    cpgsci(3); cpgsls(2);
	    cpgmove(wl[fbsp],plenv.ymin[3]); cpgdraw(wl[fbsp],plenv.ymax[3]);
	    cpgmove(wl[fbep],plenv.ymin[3]); cpgdraw(wl[fbep],plenv.ymax[3]);
	    cpgsls(1);
	    cpgmove(wl[fsp],plenv.ymin[3]); cpgdraw(wl[fsp],plenv.ymax[3]);
	    cpgmove(wl[fep],plenv.ymin[3]); cpgdraw(wl[fep],plenv.ymax[3]);
	    /* Alter fit by selecting new bounds, changing order or
	       selecting/deselecting points */
	    cpgsci(5); cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	    if (!strncmp(pgch,"?",1)) {
	      fprintf(stderr,"\
Options for windowing new plotting region:\n\
 Left mouse  : Select left- or right-limit or smoothing line to change.\n\
                 Subsequently press 'r' to register selected line to\n\
                 corresponding limit or smoothing line. Any other key to place\n\
                 line at new position.\n\
 Right mouse : Clip pixels from subsequent continuum fit.\n\
                 Right mouse again selects second corner of clip window,\n\
                 any other entry to abort selection.\n\
 Middle mouse: Unclip pixels from subsequent continuum fit.\n\
                 Middle mouse again selects second corner of clip window,\n\
                 any other entry to abort selection.\n\
 +           : Increase fit order by 1.\n\
 -           : Decrease fit order by 1.\n\
 a           : Accept fit\n\
 c           : Cancel any and all clip and unclip regions\n\
 l, u        : Alter lower/upper sigma-clipping rejection threshold.\n\
                 Points more than %5.2lf-sigma below & %5.2lf-sigma above curve\n\
                 are rejected from next iteration.\n\
 t           : Alter fit type. 1=polynomial, 2=Chebyshev, 3=Legendre.\n\
                 Current fit type is %d.\n\
 o           : Alter fit order. Current fit order is %d.\n\
 q           : Quit fitting procedure\n\n",
		      rejsigl,rejsigu,fit_typ,fit_ord);
	      cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
	    }
	    if (!strncmp(pgch,"A",1)) {
	      clickx=0.01*(plenv.xmax[1]-plenv.xmin[1]);
	      clicknp=idxfval(&(wl[sp]),cnp,wl[sp]+clickx)+2;
	      if (pgx>wl[fsp]-clickx && pgx<wl[fsp]+clickx) {
		cpgsci(15);
		cpgmove(wl[fsp],plenv.ymin[3]); cpgdraw(wl[fsp],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fsp=fbsp-1; fbsnp=1; }
		else {
		  fsp=((fsp=idxfval(&(wl[sp]),cnp,pgx))==-1) ? fbep-clicknp :
		    fsp+sp; fsp=MIN(fsp,(fbep-clicknp)); fnp=fep-fsp+1;
		    cpgsci(15); cpgmove(wl[fbsp],plenv.ymin[3]);
		    cpgdraw(wl[fbsp],plenv.ymax[3]);
		    fbsp=MAX(fbsp,(MIN((fbep-1),(fsp+fbsnp-1)))); fbsnp=fbsp-fsp+1;
		}
	      }
	      else if (pgx>wl[fep]-clickx && pgx<wl[fep]+clickx) {
		cpgsci(15);
		cpgmove(wl[fep],plenv.ymin[3]); cpgdraw(wl[fep],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fep=fbep+1; fbenp=1; }
		else {
		  fep=((fep=idxfval(&(wl[sp]),cnp,pgx))==-1) ? ep : fep+sp;
		  fep=MAX(fep,(fbsp+clicknp)); fnp=fep-fsp+1;
		  cpgsci(15); cpgmove(wl[fbep],plenv.ymin[3]);
		  cpgdraw(wl[fbep],plenv.ymax[3]);
		  fbep=MIN(fbep,(MAX((fbsp+1),(fep-fbenp+1)))); fbenp=fep-fbep+1;
		}
	      }
	      else if (pgx>wl[fbsp]-clickx && pgx<wl[fbsp]+clickx) {
		cpgsci(15); cpgmove(wl[fbsp],plenv.ymin[3]);
		cpgdraw(wl[fbsp],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fbsp=fsp+1; fbsnp=1; }
		else {
		  fbsp=((fbsp=idxfval(&(wl[sp]),cnp,pgx))==-1) ? fbep-1: fbsp+sp;
		  fbsp=MAX((fsp+1),(MIN(fbsp,(fbep-1)))); fbsnp=fbsp-fsp+1;
		}
	      }
	      else if (pgx>wl[fbep]-clickx && pgx<wl[fbep]+clickx) {
		cpgsci(15); cpgmove(wl[fbep],plenv.ymin[3]);
		cpgdraw(wl[fbep],plenv.ymax[3]);
		cpgsci(3); cpgband(6,0,0.0,0.0,&pgx,&pgy,pgch);
		if (!strncmp(pgch,"r",1)) { fbep=fep-1; fbenp=1; }
		else {
		  fbep=((fbep=idxfval(&(wl[sp]),cnp,pgx))==-1) ? fep-1 : fbep+sp;
		  fbep=MIN((fep-1),(MAX(fbep,(fbsp+1)))); fbenp=fep-fbep+1;
		}
	      }
	    }
	    else if (!strncmp(pgch,"X",1)) {
	      /* Reject pixels from next fit */
	      cpgsci(2); pgxo=pgx; pgyo=pgy;
	      cpgband(2,0,pgx,pgy,&pgx,&pgy,pgch);
	      if (!strncmp(pgch,"X",1)) {
		if (pgxo>pgx) { FSWAP(pgxo,pgx); }
		if (pgyo>pgy) { FSWAP(pgyo,pgy); }
		for (i=0,j=sp,selpix=0; i<cnp; i++,j++) {
		  if (fst[i]==1 && (wl[j]>pgxo && wl[j]<pgx) &&
		      (fl[j]>pgyo && fl[j]<pgy)) {
		    fst[i]=PCLIP; selpix++;
		  }
		}
		/* Allocate memory for a new pair of xy pairs */
		if (!nxyp && selpix) {
		  if (!(xyp=(twoxyp *)malloc((size_t)(sizeof(twoxyp)))))
		    errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\txyp array of size 1");
		}
		else if (selpix) {
		  if (!(xyp=(twoxyp *)
			realloc(xyp,(size_t)((nxyp+1)*sizeof(twoxyp)))))
		    errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\txyp array to size %d",nxyp+1);
		}
		if (selpix) {
		  xyp[nxyp].x1=pgxo; xyp[nxyp].x2=pgx;
		  xyp[nxyp].y1=pgyo; xyp[nxyp].y2=pgy;
		  xyp[nxyp++].i=0;
		}
	      }
	    }
	    else if (!strncmp(pgch,"D",1)) {
	      /* Include pixels in next fit */
	      cpgsci(2); pgxo=pgx; pgyo=pgy;
	      cpgband(2,0,pgx,pgy,&pgx,&pgy,pgch);
	      if (!strncmp(pgch,"D",1)) {
		if (pgxo>pgx) { FSWAP(pgxo,pgx); }
		if (pgyo>pgy) { FSWAP(pgyo,pgy); }
		for (i=0,j=sp,selpix=0; i<cnp; i++,j++) {
		  if (fst[i]==PCLIP && (wl[j]>pgxo && wl[j]<pgx) &&
		      (fl[j]>pgyo && fl[j]<pgy)) {
		    fst[i]=1; selpix++;
		  }
		}
		/* Allocate memory for a new pair of xy pairs */
		if (!nxyp && selpix) {
		  if (!(xyp=(twoxyp *)malloc((size_t)(sizeof(twoxyp)))))
		    errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\txyp array of size 1");
		}
		else if (selpix) {
		  if (!(xyp=(twoxyp *)
			realloc(xyp,(size_t)((nxyp+1)*sizeof(twoxyp)))))
		    errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\txyp array to size %d",nxyp+1);
		}
		if (selpix) {
		  xyp[nxyp].x1=pgxo; xyp[nxyp].x2=pgx;
		  xyp[nxyp].y1=pgyo; xyp[nxyp].y2=pgy;
		  xyp[nxyp++].i=1;
		}
	      }
	    }
	    else if (!strncmp(pgch,"c",1)) {
	      for (i=0,j=sp; i<cnp; i++,j++) fst[i]=cspec->st[j];
	      if (nxyp) { free(xyp); nxyp=0; }
	    }
	    else if (!strncmp(pgch,"t",1)) {
	      firstfit=1; while (firstfit || fit_typ<FITPOL || fit_typ>FITLEG) {
		get_input("Continuum fit type? [1=Poly; 2=Cheb; 3=Lege]","%d",
			  &fit_typ);
		if (fit_typ<FITPOL || fit_typ>FITLEG)
		  warnmsg("UVES_plot_cspec(): Fit type must be >=%d and <=%d",
			  FITPOL,FITLEG);
		firstfit=0;
	      }
	    }
	    else if (!strncmp(pgch,"o",1)) {
	      firstfit=1; while (firstfit || (fit_ord<1 || fnp<=fit_ord)) {
		get_input("Order of continuum fit?","%d",&fit_ord);
		if (fit_ord<=0)
		  warnmsg("UVES_plot_cspec(): Fit order must be > 0");
		if (fnp<=fit_ord)
		  warnmsg("UVES_plot_cspec(): Not enough valid data points\n\
\t(=%d) in selected area to fit an %d-order polynomial",fnp,fit_ord);
		firstfit=0;
	      }
	    }
	    else if (!strncmp(pgch,"+",1)) {
	      fit_ord++;
	      if (fnp<=fit_ord) {
		warnmsg("UVES_plot_cspec(): Not enough valid data points\n\
\t(=%d) in selected area to fit an %d-order polynomial",fnp,fit_ord);
		fit_ord--;
	      }
	    }
	    else if (!strncmp(pgch,"-",1)) {
	      fit_ord--;
	      if (fit_ord<1) {
		warnmsg("UVES_plot_cspec(): Fit order must be > 0");
		fit_ord++;
	      }
	    }
	    else if (!strncmp(pgch,"u",1)) {
	      firstfit=1; while (firstfit || rejsigu<DRNDTOL) {
		get_input("Upper rejection threshold in sigma?","%lf",&rejsigu);
		if (rejsigu<DRNDTOL)
		  warnmsg("UVES_plot_cspec(): Rej.-threshold must be > %5.2lg",
			  DRNDTOL);
		firstfit=0;
	      }
	    }
	    else if (!strncmp(pgch,"l",1)) {
	      firstfit=1; while (firstfit || rejsigl<DRNDTOL) {
		get_input("Lower rejection threshold in sigma?","%lf",&rejsigl);
		if (rejsigl<DRNDTOL)
		  warnmsg("UVES_plot_cspec(): Rej.-threshold must be > %5.2lg",
			  DRNDTOL);
		firstfit=0;
	      }
	    }
	  }
	  if (!strncmp(pgch,"a",1) ||
	      (par->replay && rp->exit==3 && !strncmp(pgch,"q",1))) {
	    /* Remember action by filling in an action report */
	    if (par->replay && rp->exit==3) {
	      (*act)[rp->cact].xyp=xyp; xyp=NULL;
	      (*act)[rp->cact].nxyp=nxyp; nxyp=0;
	      (*act)[rp->cact].d[0]=cspec->wl[fsp];
	      (*act)[rp->cact].d[1]=cspec->wl[fep];
	      (*act)[rp->cact].d[2]=cspec->wl[fbsp];
	      (*act)[rp->cact].d[3]=cspec->wl[fbep];
	      (*act)[rp->cact].d[4]=rejsigl; (*act)[rp->cact].d[5]=rejsigu;
	      (*act)[rp->cact].i[0]=fit_typ; (*act)[rp->cact].i[1]=fit_ord;
	    } else {
	      if (!(*nact)) {
		if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
		  errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
	      } else if (!(*act=(action *)
			   realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
		errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
	      /* If in paused replay mode then shift later actions to end */
	      if (par->replay && rp->exit==2) {
		for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
		cact=rp->cact++; (*nact)++;
	      } else cact=(*nact)++;
	      (*act)[cact].act=FCACT; (*act)[cact].xyp=xyp; xyp=NULL;
	      (*act)[cact].nxyp=nxyp; nxyp=0; (*act)[cact].d[0]=cspec->wl[fsp];
	      (*act)[cact].d[1]=cspec->wl[fep]; (*act)[cact].d[2]=cspec->wl[fbsp];
	      (*act)[cact].d[3]=cspec->wl[fbep]; (*act)[cact].d[4]=rejsigl;
	      (*act)[cact].d[5]=rejsigu; (*act)[cact].i[0]=fit_typ;
	      (*act)[cact].i[1]=fit_ord; (*act)[cact].rcmb=0;
	      (*act)[cact].nordact=1; (*act)[cact].val=1;
	      /* Alter relevant pixel values */
	      if (par->thar<=1) {
		for (i=0,j=sp; i<cnp; i++,j++) {
		  cspec->oco[j]=cspec->co[j]; scale=cspec->co[j]/fco[i];
		  fscale=(float)scale; cspec->co[j]=fco[i]; cspec->no[j]*=scale;
		  cspec->ne[j]*=scale; cspec->nf[j]*=scale; co[j]=fco[i];
		  no[j]*=fscale; ne[j]*=fscale;
		  for (k=0; k<ncon; k++) {
		    if (con[k][2]<=j && con[k][3]>=j) {
		      idx=j-spec[con[k][0]].or[con[k][1]].csidx;
		      spec[con[k][0]].or[con[k][1]].rdco[idx]=cspec->co[j];
		    }
		  }
		}
	      }
	      else {
		for (i=0,j=sp; i<cnp; i++,j++) {
		  cspec->oco[j]=cspec->co[j]; scale=cspec->co[j]-fco[i];
		  fscale=(float)scale; cspec->co[j]=fco[i]; cspec->no[j]+=scale;
		  co[j]=fco[i]; no[j]+=fscale;
		  for (k=0; k<ncon; k++) {
		    if (con[k][2]<=j && con[k][3]>=j) {
		      idx=j-spec[con[k][0]].or[con[k][1]].csidx;
		      spec[con[k][0]].or[con[k][1]].rdco[idx]=cspec->co[j];
		    }
		  }
		}
	      }
	    }
	  }
	  con_redo=0;
	  /* Clean up temporary arrays */
	  free(fco); free(fst);
	  if (par->replay && rp->exit==3) { rp->exit=0; break; }
	  if (!strncmp(pgch,"q",1)) {
	    if (nxyp) { free(xyp); nxyp=0; }
	    for (i=fsp; i<=fep; i++) co[i]=cspec->co[i];
	  }
	}
      }
    } else if (!strncmp(pgch,"N",1)) {
      /* Create a new continuumm for entire spectrum based on
	 command-line parameters, including combmeth */
      sprintf(query,"y");
      get_input("Are you sure you want to re-fit the continuum?","%s",query);
      if (strncmp(query,"n",1)) {
	if (!par->combmeth) {
	  /* Save old continuum for undo information */
	  for (i=0; i<np; i++) cspec->oco[i]=cspec->co[i];
	  /* Make new continuum */
	  if (!UVES_cspec_cont(spec,nspec,cspec,par)) {
	    nferrormsg("UVES_plot_cspec(): Error returned from\n\
\tUVES_cspec_cont() when attempting to reset continuum"); return 0;
	  }
	} else {
	  /* Reset status array before rederiving it */
	  for (i=0; i<nspec; i++) {
	    for (j=0; j<spec[i].nor; j++) {
	      if (spec[i].or[j].nuse>=MINUSE) {
		for (k=0; k<spec[i].or[j].nrdp; k++)
		  if (spec[i].or[j].rdst[k]==LCLIP && spec[i].or[j].rdst[k]==SCLIP)
		    spec[i].or[j].rdst[k]=1;
	      }
	    }
	  }
	  /* Make a rough continuum for each order of every spectrum */
	  if (!UVES_order_cont(spec,nspec,par))
	    errormsg("UVES_plot_cspec(): Unknown error returned from\n\
\tUVES_order_cont()");
	  /* Make a unified continuum for the combined spectrum */
	  if (!UVES_combine_cont(spec,nspec,cspec,1,par))
	    errormsg("UVES_plot_cspec(): Unknown error returned from\n\
\tUVES_combine_cont()");
	  /* Set flags for recombination of data */
	  cp->rccsp=1;
	  /* Set wavelength and plotting limits in combined plotting default
	     structure */
	  cp->sp=sp; cp->ep=ep; cp->swl=cspec->wl[sp]; cp->ewl=cspec->wl[ep];
	  cp->np=cnp; cp->ymn=plenv.ymin[3]; cp->ymx=plenv.ymax[3];
	  /* Set plotting modes */
	  cp->plval=plval; cp->exit=1;
	}
	/* Remember action by filling in an action report */
	if (!par->replay || rp->exit!=3) {
	  if (!(*nact)) {
	    if (!(*act=(action *)malloc((size_t)(sizeof(action)))))
	      errormsg("UVES_plot_cspec(): Cannot allocate memory for\n\
\tact array of size 1");
	  } else if (!(*act=(action *)
		       realloc(*act,(size_t)((*nact+1)*sizeof(action)))))
	    errormsg("UVES_plot_cspec(): Cannot increase memory for\n\
\tact array to size %d",*nact+1);
	  /* If in paused replay mode then shift later actions to end */
	  if (par->replay && rp->exit==2) {
	    for (j=*nact-1; j>=rp->cact; j--) (*act)[j+1]=(*act)[j];
	    cact=rp->cact++; (*nact)++;
	  } else cact=(*nact)++;
	  (*act)[cact].act=NCACT; (*act)[cact].rcmb=0;
	  (*act)[cact].nordact=1; (*act)[cact].val=1;
	}
	if (par->combmeth) (*act)[cact].rcmb=1;
	if (par->replay && rp->exit==3) { rp->exit=0; break; }
	/* Refill continuum plotting array */
	/* Fill plotting arrays */
	for (i=0; i<np; i++) {
	  co[i]=cspec->co[i]; no[i]=cspec->no[i]; ne[i]=cspec->ne[i];
	}
	con_redo=0;
      }
    } else if (!strncmp(pgch,"u",1) && cact>-1) {
      /* Undo previous action */
      if (!(*act)[cact].nordact)
	warnmsg("UVES_plot_cspec(): You cannot undo actions\n\
\tprevious to the last one you performed");
      else if (!(status=UVES_undo_lastact(spec,nspec,cspec,*act,cact+1,par))) {
	nferrormsg("UVES_plot_cspec(): Error returned from\n\
\tUVES_undo_lastact()"); return 0;
      } else if (status==2)
	warnmsg("UVES_plot_cspec(): Warning issued by UVES_undo_lastact()");
      else {
	/* If in paused replay mode get rid of undone actions */
	if (par->replay && rp->exit==2) {
	  k=(*act)[cact].nordact;
	  for (j=cact-k+1; j<=cact; j++) (*act)[j]=(*act)[j+k];
	  rp->cact-=k; (*nact)-=k; cact=rp->cact; sclinit=-1;
	} else { (*nact)-=(*act)[cact].nordact; cact=*nact-1; }
	(*act)[cact].nordact=0;
	if (!plval) {
	  for (i=0; i<np; i++) {
	    co[i]=cspec->co[i]; no[i]=cspec->no[i]; ne[i]=cspec->ne[i];
	  }
	}
      }
      con_redo=0;
    }
    else if (!strncmp(pgch,"c",1)) {
      get_input("Sig-clipping level for autorescaling orders [in sigma]?","%lf",
		&(cp->scalclip));
      con_redo=0;
    }
    else if (!strncmp(pgch,"e",1)) {
      get_input("Relative error threshold for autorescaling orders [%%/100]?",
		"%lf",&(cp->scalerr));
      con_redo=0;
    }
    else if (!strncmp(pgch," ",1)) {
      /* Check to make sure correct exit code is set */
      cp->refre=1;
      /* Set wavelength and plotting limits in combined plotting default
	 structure */
      cp->sp=sp; cp->ep=ep; cp->swl=cspec->wl[sp]; cp->ewl=cspec->wl[ep];
      cp->np=cnp; cp->ymn=plenv.ymin[3]; cp->ymx=plenv.ymax[3];
      /* Set plotting modes */
      cp->plval=plval; cp->exit=1;
    }
    else if (!strncmp(pgch,"a",1) && !par->combmeth) {
      /* Check to make sure correct exit code is set */
      if (cp->rscsp) {
	/* Set wavelength and plotting limits in combined plotting default
	   structure */
	cp->sp=sp; cp->ep=ep; cp->swl=cspec->wl[sp]; cp->ewl=cspec->wl[ep];
	cp->np=cnp; cp->ymn=plenv.ymin[3]; cp->ymx=plenv.ymax[3];
	/* Set plotting modes */
	cp->plval=plval; cp->exit=1;
      }
    }
    else if (!strncmp(pgch,"d",1)) {
      /* Allow user to delete or replace spectra in combination */
      if (!UVES_select_subspec(spec,nspec,SPSIZE*plenv.wwidth,plenv.wasp/SPSIZE,
			       par))
	errormsg("UVES_plot_cspec(): Error returned from UVES_plot_cspec()");
      cpgslct(cp->pgid);
      /* Check to make sure correct exit code is set */
      cp->rccsp=1;
      /* Set wavelength and plotting limits in combined plotting default
	 structure */
      cp->sp=sp; cp->ep=ep; cp->swl=cspec->wl[sp]; cp->ewl=cspec->wl[ep];
      cp->np=cnp; cp->ymn=plenv.ymin[3]; cp->ymx=plenv.ymax[3];
      /* Set plotting modes */
      cp->plval=plval; cp->exit=1;
    }
    else if (!strncmp(pgch,"r",1)) {
      /* Check to make sure correct exit code is set */
      if (cp->rccsp) {
	cp->rscsp=0;
	/* Set recombination code for current action */
	if (*nact) (*act)[cact].rcmb=1;
	/* Set wavelength and plotting limits in combined plotting default
	   structure */
	cp->sp=sp; cp->ep=ep; cp->swl=cspec->wl[sp]; cp->ewl=cspec->wl[ep];
	cp->np=cnp; cp->ymn=plenv.ymin[3]; cp->ymx=plenv.ymax[3];
	/* Set plotting modes */
	cp->plval=plval; cp->exit=1;
      }
    }
    else if (!strncmp(pgch,"q",1)) {
      cp->rccsp=cp->rscsp=0;
      if (par->replay && rp->exit==2) {
	/* User has paused in action replay so control should go back
	   to the replay cotroller */
	rp->actinit=rp->exit=0; break;
      } else cp->exit=1;
    }
    else con_redo=0;
  }

  /* Clean up */
  free(shwl); free(shfl); free(or); free(*con); free(con); free(*scon); free(scon);
  free(rank); free(wl); free(fl); free(er); free(co); free(no); free(ne);
  free(csq); free(ccsq); free(ncb); free(nccb);

  return 1;

}
