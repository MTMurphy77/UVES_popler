/****************************************************************************
* Initialize plotting environment
****************************************************************************/

#include <string.h>
#include <stdio.h>
#include "UVES_popler.h"

int UVES_pgenv_init(plotenv *plenv, cplot *cp) {

  float   vpy;
  
  plenv->wwidth=W_WIDTH;
  plenv->wasp=W_ASP;
  plenv->ch=C_H;
  plenv->lw=L_W;
  plenv->nxsub=N_X_SUB;
  plenv->nysub=N_Y_SUB;
  plenv->vpu=VPU;
  plenv->vpd=VPD;
  plenv->vpl=VPL;
  plenv->vpr=VPR;

  /* Arrow head style */
  cpgsah(1,45.0,0.3);

  /** Combined plot defaults **/
  if (cp!=NULL) {
    /* Set up plot geometry */
    vpy=plenv->vpu-plenv->vpd;
    cp->vpd0=plenv->vpd+VPD0*vpy; cp->vpu0=plenv->vpd+VPU0*vpy;
    cp->vpd1=plenv->vpd+VPD1*vpy; cp->vpu1=plenv->vpd+VPU1*vpy;
    cp->vpd2=plenv->vpd+VPD2*vpy; cp->vpu2=plenv->vpd+VPU2*vpy;
    cp->vpd3=plenv->vpd+VPD3*vpy; cp->vpu3=plenv->vpd+VPU3*vpy;
    cp->vpd4=plenv->vpd+VPD4*vpy; cp->vpu4=plenv->vpd+VPU4*vpy;
    /* Set wavelength and plotting limits */
    cp->swl=cp->ewl=0.0; cp->ymx=cp->ymn=0.0;
    cp->sp=cp->ep=0; cp->np=CP_NP;
    /* Set degradation factor for low resolution plot */
    cp->sfac=SFAC;
    /* Set plotting modes */
    cp->irank=cp->plval=0;
    /* Initialise exit codes */
    cp->exit=0; cp->refre=cp->rccsp=cp->rscsp=0;
    /* Initialise auto-rescaling parameters */
    cp->scalclip=cp->scalerr=0.0;
    sprintf(plenv->xlab[0],"%s","Vacuum-heliocentric wavelength [\\A]");
    sprintf(plenv->ylab[1],"%s","N\\dpix\\u");
    sprintf(plenv->ylab[2],"%s","\\gx\\d\\gn\\u\\u2\\d");
    sprintf(plenv->ylab[3],"%s","Flux");
    sprintf(plenv->ylab[4],"%s","Normalized flux");
  }

  return 1;
}
