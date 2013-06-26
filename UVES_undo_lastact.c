/****************************************************************************
* Undo the last nordact actions in the action array. Note that if the
* user has performed a recombination after an action then, if that
* action in undone, the recombination flag is moved to the previous
* action.
****************************************************************************/

#include "UVES_popler.h"
#include "error.h"

int UVES_undo_lastact(spectrum *spec, int nspec, cspectrum *cspec, action *act,
		      int nact, params *par) {

  double    scale=0;
  int       op1=0,op2=0,cp1=0,cp2=0;
  int       i=0,j=0,k=0,l=0,m=0,n=0;

  /* Set j to last action */
  j=nact-1;

  /* First check if a recombination was performed after the last
     action. If so, issue warning to the user that no undo is
     possible */
  if (act[j].rcmb==1) {
    warnmsg("UVES_undo_lastact(): Cannot undo last action\n\
\tbecause you have recombined the spectra since performing it");
    return 2;
  }

  /* Make sure last action was a valid one */
  if (!act[j].val) {
    warnmsg("UVES_undo_lastact(): Cannot undo last action\n\
\tbecause it has been rendered invalid, probably because it was skipped\n\
\tin replay mode");
    return 2;
  }

  /* Make sure that the number of actions performed in one go was
     flagged in action structure of last action */
  if (!act[j].nordact) {
    nferrormsg("UVES_undo_lastact(): Number of actions performed\n\
\t'in one go' has not been registered as variable nordact within action %d",nact);
    return 0;
  }
  
  /* Loop backwards over actions which were performed in one go */
  for (i=0,k=j; i<act[j].nordact; i++,k--) {
    /* Some actions require finding pixels on the combined spectrum ... */
    if (act[k].act<=SOACT) {
      if ((cp1=idxdval(cspec->rwl,cspec->np,act[k].d[0]))==-1)
	errormsg("UVES_undo_lastact(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",k+1,act[k].d[0]);
      if ((cp2=idxdval(&(cspec->rwl[cp1]),cspec->np-cp1,act[k].d[1]))==-1) {
	if (act[k].d[1]>cspec->rwl[cspec->np-1]) cp2=cspec->np-1-cp1;
	else
	  errormsg("UVES_undo_lastact(): Action %d: Cannot find pixel near\n\
\twavelength %lf in combined spectrum",k+1,act[k].d[1]);
      }
      cp2+=cp1;
      /* Some actions require finding pixels on individual orders as well */
      if (act[k].act==COACT || act[k].act==UOACT || act[k].act==FOACT ||
	  act[k].act==IOACT) {
	if (cp2>=spec[act[k].i[0]].or[act[k].i[1]].csidx &&
	    cp1<=spec[act[k].i[0]].or[act[k].i[1]].ceidx) {
	  cp1=MAX(cp1,spec[act[k].i[0]].or[act[k].i[1]].csidx);
	  cp2=MIN(cp2,spec[act[k].i[0]].or[act[k].i[1]].ceidx);
	}
	else
	  errormsg("UVES_undo_lastact(): Action %d: Wavelength range\n\
\tof clip-window, %lf-%lf, does not overlap with that\n\
\tof order %d in file %d,\n\t%s\n",k+1,act[k].d[0],act[k].d[1],
		   act[k].i[1]+1,act[k].i[0]+1,spec[act[k].i[0]].file);
	op1=cp1-spec[act[k].i[0]].or[act[k].i[1]].csidx;
	op2=cp2-spec[act[k].i[0]].or[act[k].i[1]].csidx;
      }
    }
    /* Enter a decision tree depending on the nature of last action */
    switch (act[k].act) {
    case COACT: 
      /** Last action clipped pixels from an order **/
    case UOACT: 
      /** Last action unclipped pixels from an order **/
      /* Replace old flags in clipped region */
      for (l=op1; l<=op2; l++) {
	if (spec[act[k].i[0]].or[act[k].i[1]].rdfl[l]>act[k].d[2] &&
	    spec[act[k].i[0]].or[act[k].i[1]].rdfl[l]<act[k].d[3])
	  spec[act[k].i[0]].or[act[k].i[1]].rdst[l]=
	    spec[act[k].i[0]].or[act[k].i[1]].ordst[l];
      }
      break;
    case CCACT:
      /** Last action clipped pixels from combined spectrum **/
      /* Replace old flags in clipped region */
      for (l=cp1; l<=cp2; l++) {
	if (cspec->fl[l]>act[k].d[2] && cspec->fl[l]<act[k].d[3]) {
	  cspec->st[l]=cspec->ost[l]; cspec->no[l]=cspec->fl[l]/cspec->co[l];
	  cspec->ne[l]=cspec->er[l]/cspec->co[l];
	  cspec->nf[l]=cspec->ef[l]/cspec->co[l];
	}
      }
      break;
    case UCACT:
      /** Last action unclipped pixels from combined spectrum **/
      /* Replace old flags in clipped region */
      for (l=cp1; l<=cp2; l++) {
	if (cspec->fl[l]>act[k].d[2] && cspec->fl[l]<act[k].d[3]) {
	  cspec->st[l]=cspec->ost[l]; cspec->no[l]=1.0;
	  cspec->ne[l]=cspec->nf[l]=-INFIN;
	}
      }
      break;
    case FOACT:
      /** Last action fitted new continuum to single order **/
    case IOACT:
      /** Last action interpolated new continuum to single order **/
      op2=(MIN(op2,(spec[act[k].i[0]].or[act[k].i[1]].nrdp-1)));
      /* Return arrays to their original values */
      if (par->thar<=1) {
	for (l=op1; l<=op2; l++) {
	  spec[act[k].i[0]].or[act[k].i[1]].rdco[l]/=
	    (scale=spec[act[k].i[0]].or[act[k].i[1]].ordco[l]);
	  spec[act[k].i[0]].or[act[k].i[1]].rdfl[l]/=scale;
	  spec[act[k].i[0]].or[act[k].i[1]].rder[l]/=scale;
	  spec[act[k].i[0]].or[act[k].i[1]].rdef[l]/=scale;
	  spec[act[k].i[0]].or[act[k].i[1]].rdme[l]/=scale;
	}
      }
      else {
	for (l=op1; l<=op2; l++)
	  spec[act[k].i[0]].or[act[k].i[1]].rdco[l]=
	    spec[act[k].i[0]].or[act[k].i[1]].ordco[l];
      }
      break;
    case FCACT:
      /** Last action fitted new continuum to combined spectrum **/
    case ICACT:
      /** Last action interpolated new continuum to combined spectrum **/
      /* Return continuum to its original value */
      for (l=cp1; l<=cp2; l++) {
	if (par->thar<=1 && cspec->co[l]!=0.0) {
	  scale=cspec->oco[l]/cspec->co[l]; cspec->co[l]=cspec->oco[l];
	  cspec->no[l]/=scale; cspec->ne[l]/=scale; cspec->nf[l]/=scale;
	} else if (par->thar==2) {
	  scale=cspec->oco[l]-cspec->co[l]; cspec->co[l]=cspec->oco[l];
	  cspec->no[l]-=scale;
	}
	/* Find contributing pixels to this combined pixel */
	for (m=0; m<nspec; m++) {
	  for (n=0; n<spec[m].nor; n++) {
	    if (spec[m].or[n].nuse>=MINUSE && l>=spec[m].or[n].csidx &&
		l<=spec[m].or[n].ceidx)
	      spec[m].or[n].rdco[l-spec[m].or[n].csidx]=cspec->co[l];
	  }
	}
      }
      break;
    case SOACT:
      /** Last action scaled a single order **/
      /* Scale arrays back to their original values */
      for (l=0; l<spec[act[k].i[0]].or[act[k].i[1]].nrdp; l++) {
	spec[act[k].i[0]].or[act[k].i[1]].rdfl[l]/=act[k].d[0];
	spec[act[k].i[0]].or[act[k].i[1]].rder[l]/=act[k].d[0];
	spec[act[k].i[0]].or[act[k].i[1]].rdef[l]/=act[k].d[0];
	spec[act[k].i[0]].or[act[k].i[1]].rdme[l]/=act[k].d[0];
      }
      break;
    case NCACT:
      /** Last action replaced continuum with new estimate **/
      /* Replace continuum with old continuum */
      for (l=0; l<cspec->np; l++) cspec->co[l]=cspec->oco[l];
      break;
    }
  }

  return 1;

}
