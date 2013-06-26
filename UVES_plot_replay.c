/****************************************************************************
* Plot replay control window
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "UVES_popler.h"
#include "input.h"
#include "error.h"

int UVES_plot_replay(rplot *rp, action *act, int nact, plotbut *but, int nbut,
		     params *par) {

  int     i=0;
  char    key1[LNGSTRLEN]="\0",key2[LNGSTRLEN]="\0";
  plotenv plenv;
  extern  char   *progname;

  /** Initialize plotting **/
  if (!rp->init) {
    /* If plot has not been initialized then initialize it */
    UVES_pgenv_init(&plenv,NULL); plenv.wwidth*=RPSIZE; plenv.wasp=1.5;
    plenv.vpl=plenv.vpd=0.00; plenv.vpr=plenv.vpu=1.00;
    sprintf(rp->ptitle,"%s: Replay control panel",progname); cpgask(0);
    pg_open(&plenv,"/XWIN",rp->ptitle,1); cpgqid(&(rp->pgid)); cpgpage();
    rp->zx=rp->zxmin=0.5; rp->zzx=5.0*rp->zx; rp->zymin=0.01; rp->zy=5.0*rp->zymin;
    rp->zzy=5.0*rp->zy;
    rp->actinit=rp->exun=rp->skne=rp->zomo=rp->decx=rp->prre=0;
    rp->init=rp->decy=1;
  } else {
    /* If plot has already been initialized then select correct PGPLOT window */
    cpgslct(rp->pgid); plenv.vpl=plenv.vpd=0.00; plenv.vpr=plenv.vpu=1.00;
  }

  /* Begin plotting buffer */
  cpgbbuf(); cpgeras();

  /* Describe action to be performed */
  switch (act[rp->cact].act) {
  case COACT:
    sprintf(key1,"Act. %d: %d: Clip pixels from",rp->cact+1,
	    act[rp->cact].act);
    sprintf(key2,"order %d of spectrum %d",act[rp->cact].i[1]+1,
	    act[rp->cact].i[0]+1);
    rp->nav=(rp->exun) ? 0 : 1; break;
  case CCACT:
    sprintf(key1,"Act. %d: %d: Clip pixels from",rp->cact+1,
	    act[rp->cact].act);
    sprintf(key2,"combined spectrum");
    rp->nav=(rp->exun) ? 0 : 1; break;
  case UOACT:
    sprintf(key1,"Act. %d: %d: Unclip pixels from",rp->cact+1,
	    act[rp->cact].act);
    sprintf(key2,"combined spectrum");
    rp->nav=(rp->exun) ? 0 : 1; break;
  case UCACT:
    sprintf(key1,"Act. %d: %d: Unclip pixels in\\n\
combined spectrum",rp->cact+1,act[rp->cact].act);
    rp->nav=(rp->exun) ? 0 : 1; break;
  case FOACT:
    sprintf(key1,"Act. %d: %d: Fit new cont.",rp->cact+1,
	    act[rp->cact].act);
    sprintf(key2,"to order %d of spectrum %d",act[rp->cact].i[1]+1,
	    act[rp->cact].i[0]+1);
    rp->nav=0; break;
  case FCACT:
    sprintf(key1,"Act. %d: %d: Fit new cont.",rp->cact+1,
	    act[rp->cact].act);
    sprintf(key2,"to combined spectrum");
    rp->nav=0; break;
  case IOACT:
    sprintf(key1,"Act. %d: %d: Interp. new cont.",rp->cact+1,
	    act[rp->cact].act);
    sprintf(key2,"to order %d of spectrum %d",act[rp->cact].i[1]+1,
	    act[rp->cact].i[0]+1);
    rp->nav=0; break;
  case ICACT:
    sprintf(key1,"Act. %d: %d: Interp. new cont.",rp->cact+1,
	    act[rp->cact].act);
    sprintf(key2,"to combined spectrum");
    rp->nav=0; break;
  case SOACT:
    sprintf(key1,"Act. %d: %d: Scale order %d of",rp->cact+1,
	    act[rp->cact].act,act[rp->cact].i[1]+1);
    sprintf(key2,"spectrum %d",act[rp->cact].i[0]+1);
    rp->nav=0; break;
  case ARACT:
    sprintf(key1,"Act. %d: %d: Rescale orders",rp->cact+1,
	    act[rp->cact].act);
    sprintf(key2,"and recombine spectrum");
    rp->nav=0; break;
  case NCACT:
    sprintf(key1,"Act. %d: %d: Create new",rp->cact+1,
	    act[rp->cact].act);
    sprintf(key2,"automatic continuum");
    rp->nav=0; break;
  }

  /* Fill in values for positions and attributes of buttons */
  sprintf(but[0].lab,"Zoom");
  but[0].x1=0.05; but[0].x2=0.38; but[0].y1=0.70; but[0].y2=0.81;
  sprintf(but[1].lab,"Move");
  but[1].x1=0.05; but[1].x2=0.38; but[1].y1=0.58; but[1].y2=0.69;
  sprintf(but[2].lab,"Accel. X");
  but[2].x1=0.62; but[2].x2=0.95; but[2].y1=0.70; but[2].y2=0.81;
  sprintf(but[3].lab,"Decel. X");
  but[3].x1=0.62; but[3].x2=0.95; but[3].y1=0.58; but[3].y2=0.69;
  sprintf(but[4].lab,"Accel. Y");
  but[4].x1=0.62; but[4].x2=0.95; but[4].y1=0.31; but[4].y2=0.42;
  sprintf(but[5].lab,"Decel. Y");
  but[5].x1=0.62; but[5].x2=0.95; but[5].y1=0.19; but[5].y2=0.30;
  sprintf(but[6].lab,"<<");
  but[6].x1=0.03; but[6].x2=0.21; but[6].y1=0.44; but[6].y2=0.56;
  sprintf(but[7].lab,"<");
  but[7].x1=0.22; but[7].x2=0.40; but[7].y1=0.44; but[7].y2=0.56;
  sprintf(but[8].lab,">");
  but[8].x1=0.60; but[8].x2=0.78; but[8].y1=0.44; but[8].y2=0.56;
  sprintf(but[9].lab,">>");
  but[9].x1=0.79; but[9].x2=0.97; but[9].y1=0.44; but[9].y2=0.56;
  sprintf(but[10].lab,"\\(0537)\\(0537)");
  but[10].x1=0.41; but[10].x2=0.59; but[10].y1=0.70; but[10].y2=0.82;
  sprintf(but[11].lab,"\\(0537)");
  but[11].x1=0.41; but[11].x2=0.59; but[11].y1=0.57; but[11].y2=0.69;
  sprintf(but[12].lab,"V");
  but[12].x1=0.41; but[12].x2=0.59; but[12].y1=0.31; but[12].y2=0.43;
  sprintf(but[13].lab,"VV");
  but[13].x1=0.41; but[13].x2=0.59; but[13].y1=0.18; but[13].y2=0.30;
  but[14].x1=0.01; but[14].x2=0.33; but[14].y1=0.09; but[14].y2=0.17;
  but[15].x1=0.34; but[15].x2=0.66; but[15].y1=0.09; but[15].y2=0.17;
  but[16].x1=0.67; but[16].x2=0.99; but[16].y1=0.09; but[16].y2=0.17;
  sprintf(but[17].lab,"Pause");
  but[17].x1=0.01; but[17].x2=0.33; but[17].y1=0.01; but[17].y2=0.08;
  sprintf(but[18].lab,"Finish");
  but[18].x1=0.34; but[18].x2=0.66; but[18].y1=0.01; but[18].y2=0.08;
  sprintf(but[19].lab,"QUIT");
  but[19].x1=0.67; but[19].x2=0.99; but[19].y1=0.01; but[19].y2=0.08;
  for (i=0; i<nbut; i++) {
    but[i].ch=1.8*plenv.ch; but[i].bgc=0; but[i].fgc=but[i].boc=5;
  }

  /* Decide on other context-dependent button names and colours */
  if (rp->nav) {
    if (rp->zomo) { but[0].fgc=but[0].boc=14; but[1].fgc=but[1].boc=5; }
    else { but[0].fgc=but[0].boc=5; but[1].fgc=but[1].boc=14; }
  } else for (i=0; i<=13; i++) but[i].fgc=but[i].boc=14;
  if (!rp->decx) but[3].fgc=but[3].boc=14;
  if (!rp->decy) but[5].fgc=but[5].boc=14;
  if (rp->exun) sprintf(but[14].lab,"Undo"); else sprintf(but[14].lab,"Execute");
  if (rp->skne) sprintf(but[15].lab,"Next"); else sprintf(but[15].lab,"Skip");

  /* Decide if it's possible to recombine the spectrum, given the
     action history */
  if (rp->prre) {
    sprintf(but[16].lab,"Recombine");
    if (rp->rcmb && act[rp->cact].act!=ARACT) {
      /* Search previous actions to see if there have been any valid
	 actions which makes recombination meaningful */
      for (i=rp->cact-1; i>=0; i--) {
	if (act[i].val && (act[i].act==ARACT || act[i].rcmb==1)) {
	  rp->prre=-1; break;
	}
	else if (act[i].val &&
		 (act[i].act==COACT || act[i].act==UOACT ||
		  act[i].act==FOACT || act[i].act==IOACT ||
		  act[i].act==SOACT)) { rp->prre=1; break; }
      }
      if (rp->prre==1) { but[16].bgc=2; but[16].fgc=3; but[16].boc=4; }
    }
    else if (act[rp->cact].act==ARACT) rp->prre=-1;
    else if (act[rp->cact].act==NCACT && par->combmeth) rp->prre=-1;
  } else sprintf(but[16].lab,"Previous");

  /* Decide whether it's possible to go back to previous action */
  but[16].ch=(rp->prre) ? 1.6*plenv.ch : 1.8*plenv.ch;
  if (!rp->cact && !rp->prre) rp->prre=-1;
  else if (rp->cact && !rp->prre && act[rp->cact-1].rcmb) rp->prre=-1;
  if (rp->prre==-1) but[16].fgc=but[16].boc=14;

  /* Plot main replay display window */
  cpgsvp(plenv.vpl,plenv.vpr,plenv.vpd,plenv.vpu); cpgswin(0.0,1.0,0.0,1.0);
  cpgsci(15); cpgrect(0.0,1.0,0.0,1.0); cpgsci(0); cpgsfs(1);
  cpgrect(0.01,0.99,0.85,0.99); cpgsci(3); cpgsfs(2);
  cpgrect(0.01,0.99,0.85,0.99); cpgsfs(1); cpgsci(7); cpgsch(1.7*plenv.ch);
  cpgtext(0.03,0.94,key1); cpgtext(0.05,0.88,key2);
  for (i=0; i<nbut; i++) {
    if (!pg_button(&(but[i]))) {
      nferrormsg("UVES_plot_replay(): Error returned from\n\
\tpg_button()"); return 0;
    }
  }

  /* Flush plotting buffer */
  cpgebuf();

  return 1;

}
