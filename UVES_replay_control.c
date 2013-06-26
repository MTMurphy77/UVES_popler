/****************************************************************************
* Bring up a new PGPLOT window which allows the user to control the
* replay of individual actions
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "UVES_popler.h"
#include "input.h"
#include "stats.h"
#include "memory.h"
#include "error.h"
#include "const.h"

#define  NBUT 20

int UVES_replay_control(spectrum *spec, int nspec, cspectrum *cspec, cplot *cp,
			rplot *rp, action **act, int *nact, int *nact_save,
			params *par) {

  double  accdec=2.0,scalfac=5.0;
  double  swl=0.0,ewl=0.0,pix=0.0,ddum=0.0;
  float   pgx=0.0,pgy=0.0;
  int     sp=0,ep=0,np=0;
  int     i=0,j=0,k=0;
  char    pgch[SHORTLEN]="p";
  statset stat;
  plotbut but[NBUT];

  /* Loop around actions */
  while (rp->cact<*nact && !rp->exit) {

    /* Select and plot the replay control window */
    if (!UVES_plot_replay(rp,*act,*nact,but,NBUT,par)) {
      nferrormsg("UVES_replay_control(): Error returned from\n\
\tUVES_plot_replay() when attempting to plot action %d",rp->cact+1); return 0;
    }

    /* Do we need to initialize the plot for the current action? */
    if (!rp->actinit) {
      /* If we are to auto-rescale the orders then we must check for
	 previous manually rescaled orders since the last
	 recombination */
      if ((*act)[rp->cact].act==ARACT) {
	for (i=rp->cact-1; i>=0; i--) {
	  if ((*act)[i].val && ((*act)[i].act==ARACT || (*act)[i].rcmb==1)) {
	    (*act)[rp->cact].val=0; rp->cact++; rp->actinit=0; rp->prre=-1; break;
	  }
	  else if ((*act)[i].val && (*act)[i].act==SOACT) break;
	}
	if (i<0) {
	  (*act)[rp->cact].val=0; rp->cact++; rp->actinit=0; rp->prre=-1;
	}
	if (rp->cact>=*nact) break;
      }
      /* If we are to form a new continuum in the order-scaling
	 combination mode then we must check to make sure a
	 recombination has been performed since we last set a new
	 continuum */
      if ((*act)[rp->cact].act==NCACT && !par->combmeth) {
	for (i=rp->cact-1; i>=0; i--) {
	  if ((*act)[i].val && (*act)[i].rcmb==1) break;
	  else if ((*act)[i].val && (*act)[i].act==NCACT) {
	    (*act)[rp->cact].val=0; rp->cact++; rp->actinit=0; rp->prre=-1; break;
	  }
	}
	if (i<0) {
	  (*act)[rp->cact].val=0; rp->cact++; rp->actinit=0; rp->prre=-1;
	}
	if (rp->cact>=*nact) break;
      }
      /* Determine the clipping area to be plotted in the Spectrum
	 Navigator for this action */
      if ((*act)[rp->cact].act<=ICACT) {
	swl=(*act)[rp->cact].d[0]-FFRAC*
	  ((*act)[rp->cact].d[1]-(*act)[rp->cact].d[0]);
	ewl=(*act)[rp->cact].d[1]+FFRAC*
	  ((*act)[rp->cact].d[1]-(*act)[rp->cact].d[0]);
	ewl=(MIN(cspec->rwl[cspec->np-1],ewl));
	if ((cp->sp=idxdval(cspec->rwl,cspec->np,swl))==-1) {
	  nferrormsg("UVES_replay_control(): Cannot find pixel with\n\
\twavelength %lf in combined spectrum",(*act)[rp->cact].d[0]); return 0;
	}
	np=cspec->np-cp->sp;
	if ((cp->ep=idxdval(&(cspec->rwl[cp->sp]),np,ewl))==-1) cp->ep=cspec->np-1;
	else (cp->ep)+=cp->sp;
      } else if ((*act)[rp->cact].act==SOACT) {
	sp=spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].csidx;
	ep=spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].ceidx;
	np=ep-sp+1;
	cp->sp=(MAX(0,(sp-(int)(FFRAC*(double)np))));
	cp->ep=(MIN((cspec->np-1),(ep+(int)(FFRAC*(double)np))));
      }
      if ((*act)[rp->cact].act<=UCACT) {
	cp->ymn=(float)((*act)[rp->cact].d[2]-FFRAC*
			((*act)[rp->cact].d[3]-(*act)[rp->cact].d[2]));
	cp->ymx=(float)((*act)[rp->cact].d[3]+FFRAC*
			((*act)[rp->cact].d[3]-(*act)[rp->cact].d[2]));
      } else if ((*act)[rp->cact].act==FOACT || (*act)[rp->cact].act==IOACT ||
		 (*act)[rp->cact].act==SOACT) {
	if ((*act)[rp->cact].act==FOACT || (*act)[rp->cact].act==IOACT) {
	  sp=MAX(spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].csidx,
		 cp->sp);
	  ep=MIN(spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].ceidx,
		 cp->ep);
	  np=ep-sp+1;
	  sp-=spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].csidx;
	} else {
	  sp=0; np=spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].np;
	}
	ep=sp+np-1;
	if (!stats(
          &(spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].rdfl[sp]),
	  &(spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].rder[sp]),
	  &(spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].rdef[sp]),
	  &(spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].rdme[sp]),
	  &(spec[(*act)[rp->cact].i[0]].or[(*act)[rp->cact].i[1]].rdst[sp]),np,2,
	  &stat)) {
	  nferrormsg("UVES_replay_control(): Error returned from stats()\n\
\twhen doing statistics between pixels %d and %d in order %d of spectrum %d\n\t%s",
		     sp+1,ep+1,(*act)[rp->cact].i[1]+1,(*act)[rp->cact].i[1]+1,
		     spec[(*act)[rp->cact].i[0]].file);
	}
	cp->ymn=(float)(stat.med-scalfac*stat.siqr);
	cp->ymx=(float)(stat.med+scalfac*stat.siqr);
	if ((*act)[rp->cact].act==FOACT) {
	  for (i=0; i<(*act)[rp->cact].nxyp; i++) {
	    cp->ymn=MIN(cp->ymn,(*act)[rp->cact].xyp[i].y1);
	    cp->ymx=MAX(cp->ymx,(*act)[rp->cact].xyp[i].y2);
	  }
	} else if ((*act)[rp->cact].act==IOACT) {
	  for (i=0; i<(*act)[rp->cact].nxyp; i++) {
	    cp->ymn=MIN(cp->ymn,(*act)[rp->cact].xyp[i].y1);
	    cp->ymx=MAX(cp->ymx,(*act)[rp->cact].xyp[i].y1);
	  }
	}
      } else if ((*act)[rp->cact].act==FCACT || (*act)[rp->cact].act==ICACT) {
	np=cp->ep-cp->sp+1;
	if (!stats(&(cspec->fl[cp->sp]),&(cspec->er[cp->sp]),&(cspec->ef[cp->sp]),
		   NULL,&(cspec->st[cp->sp]),np,2,&stat)) {
	  nferrormsg("UVES_replay_control(): Error returned from stats()\n\
\twhen doing statistics between pixels %d and %d (%lf-%lf)\n\
\tof combined spectrum",cp->sp+1,cp->ep+1,cspec->wl[cp->sp],cspec->wl[cp->ep]);
	  return 0;
	}
	cp->ymn=(float)(stat.med-scalfac*stat.siqr);
	cp->ymx=(float)(stat.med+scalfac*stat.siqr);
	if ((*act)[rp->cact].act==FCACT) {
	  for (i=0; i<(*act)[rp->cact].nxyp; i++) {
	    cp->ymn=MIN(cp->ymn,(*act)[rp->cact].xyp[i].y1);
	    cp->ymx=MAX(cp->ymx,(*act)[rp->cact].xyp[i].y2);
	  }
	} else if ((*act)[rp->cact].act==ICACT) {
	  for (i=0; i<(*act)[rp->cact].nxyp; i++) {
	    cp->ymn=MIN(cp->ymn,(*act)[rp->cact].xyp[i].y1);
	    cp->ymx=MAX(cp->ymx,(*act)[rp->cact].xyp[i].y1);
	  }
	}
      }
      cp->np=cp->ep-cp->sp+1; cp->swl=cspec->wl[cp->sp]; cp->ewl=cspec->wl[cp->ep];
      /* Rank the orders in combination rank when scaling orders */
      cp->irank=((*act)[rp->cact].act==SOACT) ? 1 : 0;
      /* Remember recombination flag of action and reset flag in action itself */
      rp->rcmb=(*act)[rp->cact].rcmb; (*act)[rp->cact].rcmb=0;
      /* Tell rest of program that this action has been initialized for plotting */
      rp->actinit=1;
      /* Decide which exit code to use. This will decide whether the
	 user will be able to navigate the spectrum and modify the
	 current action parameters on-the-fly (exit=4) or only view
	 the relevant portion of the spectrum */
      if ((*act)[rp->cact].act>=FOACT && (*act)[rp->cact].act<=ICACT &&
	  !rp->exun) rp->exit=3;
      else rp->exit=4;
      /* Exit while loop here so that the Spectrum Navigator panel can
	 be plotted, highlighting this action */
      break;
    }

    /* Allow user to make selections */
    cpgband(0,0,0.0,0.0,&pgx,&pgy,pgch);
    if (!strncmp(pgch,"?",1)) {
      fprintf(stderr,"\
Options for replay panel:\n\
 Left mouse on \"Move/Zoom\"   : Toggles between moving and zooming in/out\n\
                                 clip/unclip box\n\
 Left mouse on \"Accel/Decel\" : Increases/decreases moving/zooming increment\n\
                                 in X or Y directions\n\
 Left mouse on single arrows : Moves/zooms clip/unclip box\n\
 Left mouse on double arrows : Moves/zooms clip/unclip box by larger amount\n\
 Left mouse on \"Execute\"     : Executes current action\n\
 Left mouse on \"Undo\"        : Undoes current action\n\
 Left mouse on \"Skip\"        : Skips current action, moves to next one\n\
 Left mouse on \"Next\"        : Move to next action after executing current one\n\
 Left mouse on \"Pause\"       : Returns normal control to Spectrum Navigator\n\
                                 panel, allowing user to insert new actions\n\
 Left mouse on \"Finish\"      : Completes all remaining actions\n\
 g                           : Goto action number, executing intervening actions\n\
 s                           : Skip to action number, skipping other actions\n\
 q or left mouse on \"QUIT\"   : Quits replay WITHOUT executing remaining\n\
                                actions\n");
    }
    else if (!strncmp(pgch,"q",1)) rp->exit=1;
    else if (!strncmp(pgch,"g",1) && rp->cact<*nact-1) {
      /* Goto to order number */
      i=rp->cact+2; get_input("Goto to action number?","%d",&i); i--;
      if (i<rp->cact) warnmsg("UVES_replay_control(): Action %d invalid.\n\
\tMust go to action later than current one (=%d)",i+1,rp->cact+1);
      else if (i>=*nact) warnmsg("UVES_replay_control(): Action %d invalid.\n\
\tMust go to action earlier than final one (=%d)",i+1,*nact);
      else {
	if (rp->exun) j=i-rp->cact-1; else j=i-rp->cact;
	if (!UVES_past_actions(spec,nspec,cspec,(*act),i,-1*j,par)) {
	  nferrormsg("UVES_replay_control(): Error returned from\n\
\tUVES_past_actions() when attempting actions %d to %d",rp->cact+1,i+1); return 0;
	}
	rp->cact=i; rp->exun=rp->skne=rp->prre=rp->actinit=0;
      }
    }
    else if (!strncmp(pgch,"s",1) && rp->cact<*nact-1) {
      /* Skip to order number */
      i=rp->cact+2; get_input("Skip to action number?","%d",&i); i--;
      if (i<rp->cact) warnmsg("UVES_replay_control(): Action %d invalid.\n\
\tMust skip to action later than current one (=%d)",i+1,rp->cact+1);
      else if (i>=*nact) warnmsg("UVES_replay_control(): Action %d invalid.\n\
\tMust skip to action earlier than final one (=%d)",i+1,*nact);
      else {
	for (j=rp->cact; j<i; j++) (*act)[j].val=0;
	rp->cact=i; rp->exun=rp->skne=rp->prre=rp->actinit=0;
      }
    }
    else if (!strncmp(pgch,"A",1)) {
      /* Figure out which button was pressed, if any */
      for (i=0,j=-1; i<NBUT; i++) {
	if (pgx>=but[i].x1 && pgx<=but[i].x2 && pgy>=but[i].y1 && pgy<=but[i].y2) {
	  j=i; break;
	}
      }
      /* Make sure selected button makes sense in this context */
      k=1; if (!rp->nav && j<=13) k=0; if (rp->prre==-1 && j==16) k=0;
      if (k) {
	/* Decide on what to do given button pressed */
	if (j>=6 && j<=9) {
	  /* Work out typical pixel size in Angstroems for clip/unclip region */
	  if (par->linear) pix=par->disp;
	  else
	    pix=par->disp*0.5*((*act)[*nact].d[0]+(*act)[rp->cact].d[1])/C_C_K;
	} else if (j>=10 && j<=13)
	  ddum=(*act)[rp->cact].d[3]-(*act)[rp->cact].d[2];
	switch (j) {
	case 0:
	  /* Toggle to zoom mode */
	  rp->zomo=0;
	  break;
	case 1:
	  /* Toggle to move mode */
	  rp->zomo=1;
	  break;
	case 2:
	  /* Accelerate X move/zoom factor */
	  rp->zx*=accdec; rp->zzx=5.0*rp->zx; rp->decx=1;
	  break;
	case 3:
	  /* Decelerate X move/zoom factor */
	  if ((rp->zx/=accdec)<=rp->zxmin) { rp->zx=rp->zxmin; rp->decx=0; }
	  rp->zzx=5.0*rp->zx;
	  break;
	case 4:
	  /* Accelerate Y move/zoom factor */
	  rp->zy*=accdec; rp->zzy=5.0*rp->zy; rp->decy=1;
	  break;
	case 5:
	  /* Decelerate Y move/zoom factor */
	  if ((rp->zy/=accdec)<=rp->zymin) { rp->zy=rp->zymin; rp->decy=0; }
	  rp->zzy=5.0*rp->zy;
	  break;
	case 6:
	  /* Double move left or double zoom in */
	  (*act)[rp->cact].d[1]-=pix*rp->zzx;
	  if (rp->zomo) (*act)[rp->cact].d[0]-=pix*rp->zzx;
	  else (*act)[rp->cact].d[0]+=pix*rp->zzx;
	  rp->exit=4;
	  break;
	case 7:
	  /* Single move left or single zoom in */
	  (*act)[rp->cact].d[1]-=pix*rp->zx;
	  if (rp->zomo) (*act)[rp->cact].d[0]-=pix*rp->zx;
	  else (*act)[rp->cact].d[0]+=pix*rp->zx;
	  rp->exit=4;
	  break;
	case 8:
	  /* Single move right or single zoom out */
	  (*act)[rp->cact].d[1]+=pix*rp->zx;
	  if (rp->zomo) (*act)[rp->cact].d[0]+=pix*rp->zx;
	  else (*act)[rp->cact].d[0]-=pix*rp->zx;
	  rp->exit=4;
	  break;
	case 9:
	  /* Double move right or double zoom out */
	  (*act)[rp->cact].d[1]+=pix*rp->zzx;
	  if (rp->zomo) (*act)[rp->cact].d[0]+=pix*rp->zzx;
	  else (*act)[rp->cact].d[0]-=pix*rp->zzx;
	  rp->exit=4;
	  break;
	case 10:
	  /* Double move up or double zoom out */
	  (*act)[rp->cact].d[3]+=rp->zzy*ddum;
	  if (rp->zomo) (*act)[rp->cact].d[2]+=rp->zzy*ddum;
	  else (*act)[rp->cact].d[2]-=rp->zzy*ddum; 
	  rp->exit=4;
	  break;
	case 11:
	  /* Single move up or single zoom out */
	  (*act)[rp->cact].d[3]+=rp->zy*ddum;
	  if (rp->zomo) (*act)[rp->cact].d[2]+=rp->zy*ddum;
	  else (*act)[rp->cact].d[2]-=rp->zy*ddum; 
	  rp->exit=4;
	  break;
	case 12:
	  /* Single move down or single zoom in */
	  (*act)[rp->cact].d[3]-=rp->zy*ddum;
	  if (rp->zomo) {
	    (*act)[rp->cact].d[2]-=rp->zy*ddum; (*act)[rp->cact].d[2]+=1.e-6;
	  }
	  else (*act)[rp->cact].d[2]+=rp->zy*ddum; 
	  rp->exit=4;
	  break;
	case 13:
	  /* Double move down or double zoom in */
	  (*act)[rp->cact].d[3]-=rp->zzy*ddum;
	  if (rp->zomo) {
	    (*act)[rp->cact].d[2]-=rp->zzy*ddum; (*act)[rp->cact].d[2]+=1.e-6;
	  }
	  else (*act)[rp->cact].d[2]+=rp->zzy*ddum; 
	  rp->exit=4;
	  break;
	case 14:
	  /* Execute or undo current action */
	  if (rp->exun) {
	    /* Undo current action */
	    if (!UVES_undo_lastact(spec,nspec,cspec,&((*act)[rp->cact]),1,par)) {
	      nferrormsg("UVES_replay_control(): Error returned from\n\
\tUVES_undo_lastact() when attempting to undo action %d",rp->cact+1); return 0;
	    }
	    rp->exun=rp->skne=0; rp->prre=-1;
	    if ((*act)[rp->cact].act>=FOACT && (*act)[rp->cact].act<=ICACT)
	      rp->exit=3;
	    else rp->exit=4;
	  } else {
	    /* Execute current action */
	    (*act)[rp->cact].val=1;
	    if (!UVES_past_actions(spec,nspec,cspec,*act,rp->cact+1,-1,par)) {
	      nferrormsg("UVES_replay_control(): Error returned from\n\
\tUVES_past_actions() when attempting action %d",rp->cact+1); return 0;
	    }
	    if ((*act)[rp->cact].act==ARACT) (*act)[rp->cact].rcmb=1;
	    if ((*act)[rp->cact].act==ARACT ||
		((*act)[rp->cact].act==NCACT && par->combmeth)) {
	      rp->cact++; rp->exun=rp->skne=rp->actinit=0; rp->prre=-1;
	    }
	    else rp->exun=rp->skne=rp->prre=1;
	    if ((*act)[rp->cact].act>=FOACT || (*act)[rp->cact].act<=ICACT)
	      if (!UVES_plot_replay(rp,*act,*nact,but,NBUT,par)) {
		nferrormsg("UVES_replay_control(): Error returned from\n\
\tUVES_plot_replay() when attempting to plot action %d",rp->cact+1); return 0;
	      }
	    rp->exit=4;
	  }
	  break;
	case 15:
	  /* Skip action or move to next one */
	  /* If skipping action, remove it from the action array by
	     setting the "val" switch to zero */
	  if (!rp->skne) (*act)[rp->cact].val=0;
	  rp->cact++; rp->exun=rp->skne=rp->prre=rp->actinit=0;
	  break;
	case 16:
	  /* Go back to previous action or recombine spectrum */
	  if (rp->prre) {
	    /* Recombine all orders and spectra */
	    if (!UVES_combine_spec(spec,nspec,cspec,par))
	      errormsg("UVES_replay_control(): Error returned\n\
\tfrom UVES_combine_spec()");
	    (*act)[rp->cact].rcmb=1; rp->cact++;
	    rp->exun=rp->skne=rp->actinit=0; rp->prre=-1;
	  } else {
	    /* Go back to previous action */
	    if (rp->cact) rp->cact--;
	    if ((*act)[rp->cact].val) {
	      rp->exun=rp->skne=rp->prre=1; rp->actinit=0;
	    }
	    else rp->exun=rp->skne=rp->prre=rp->actinit=0;
	  }
	  break;
	case 17:
	  /* Pause the replay control and pass control back to
	     Spectrum Navigator */
	  (*act)[rp->cact].rcmb=rp->rcmb;
	  rp->exit=2;
	  break;
	case 18:
	  /* Finish off the rest of the actions non-interactively,
	     close replay controller and pass control back to Spectrum
	     Navigator */
	  /* Execute current action and rest of actions */
	  if (!UVES_past_actions(spec,nspec,cspec,(*act),*nact,
	      -1*(*nact-rp->cact),par)) {
	    nferrormsg("UVES_replay_control(): Error returned from\n\
\tUVES_past_actions() when attempting actions from %d onwards",rp->cact+1);
	    return 0;
	  }
	case 19: 
	  /* Quit */
	  rp->exit=1;
	  break;	  
	}
      }
    }
  }

  /* Close window if no longer needed */
  if (rp->cact==*nact) rp->exit=1;
  if (rp->exit==1) { cpgslct(rp->pgid); cpgclos(); }

  return 1;

}
