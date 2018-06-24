/****************************************************************************
* Initialize (zero) parameters set
****************************************************************************/

#include <string.h>
#include "UVES_popler.h"

int UVES_params_init(params *par) {

  par->clipsig=par->contpctl=par->contsigl=par->contsigu=par->disp=-1.0;
  par->lrgerr=par->scalclip=par->scalerr=par->ordsig=par->ordsignbr=-1.0;
  par->ordsigzero=par->ordmedfrac=par->ordmedrej=par->pctllya=par->pctlred=-1.0;
  par->rsiglyal=par->rsiglyau=par->rsigredl=par->rsigredu=par->vclya=par->vcred=-1.0;
  par->vlya=par->zem=par->version=-1.0;
  par->atmask=par->combmeth=par->contftyp=par->contord=par->contwgt=par->cordlya=-1;
  par->cordred=par->dat=par->distort=par->ftyplya=par->ftypred=par->helio=par->linear=-1;
  par->macmap=par->nocont=par->nordclip=par->nordsig=par->rankspec=par->raw=par->replay=-1;
  par->save=par->scale=par->scalmeth=par->nscalclip=par->thar=par->vacwl=par->vshift=-1;
  par->backvers=-1;
  par->filetype=-2;
  strcpy(par->prefix,"\0"); strcpy(par->macmapfile,"\0"); strcpy(par->vshiftfile,"\0");
  strcpy(par->scalefile,"\0");

  return 1;

}
