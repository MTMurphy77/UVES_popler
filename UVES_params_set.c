/****************************************************************************
* Set unset parameters to default values
****************************************************************************/

#include <stdio.h>
#include "UVES_popler.h"

int UVES_params_set(params *par) {

  if (par->clipsig==-1.0) par->clipsig=CLIPSIG;
  if (par->contpctl==-1.0) par->contpctl=CONTPCTL;
  if (par->contsigl==-1.0) par->contsigl=CONTSIGL;
  if (par->contsigu==-1.0) par->contsigu=CONTSIGU;
  if (par->disp==-1.0) par->disp=0.0;
  if (par->lrgerr==-1.0) par->lrgerr=LRGERR;
  if (par->ordsig==-1.0) par->ordsig=ORDSIG;
  if (par->ordsignbr==-1.0) par->ordsignbr=ORDSIGNBR;
  if (par->ordsigzero==-1.0) par->ordsigzero=ORDSIGZERO;
  if (par->ordmedfrac==-1.0) par->ordmedfrac=ORDMEDFRAC;
  if (par->ordmedrej==-1.0) par->ordmedrej=ORDMEDREJ;
  if (par->pctllya==-1.0) par->pctllya=PCTLLYA;
  if (par->pctlred==-1.0) par->pctlred=PCTLRED;
  if (par->rsiglyal==-1.0) par->rsiglyal=RSIGLYAL;
  if (par->rsiglyau==-1.0) par->rsiglyau=RSIGLYAU;
  if (par->rsigredl==-1.0) par->rsigredl=RSIGREDL;
  if (par->rsigredu==-1.0) par->rsigredu=RSIGREDU;
  if (par->scalclip==-1.0) par->scalclip=SCALCLIP;
  if (par->scalerr==-1.0) par->scalerr=SCALERR;
  if (par->vclya==-1.0) par->vclya=VCLYA;
  if (par->vcred==-1.0) par->vcred=VCRED;
  if (par->vlya==-1.0) par->vlya=VLYA;
  if (par->zem==-1.0) par->zem=ZEM;
  if (par->version==-1.0) par->version=VERSION;

  if (par->atmask==-1) par->atmask=0;
  if (par->combmeth==-1) par->combmeth=COMBMETH;
  if (par->contftyp==-1) par->contftyp=CONTFTYP;
  if (par->contord==-1) par->contord=CONTORD;
  if (par->contwgt==-1) par->contwgt=CONTWGT;
  if (par->cordlya==-1) par->cordlya=CORDLYA;
  if (par->cordred==-1) par->cordred=CORDRED;
  if (par->dat==-1) par->dat=DAT;
  if (par->distort==-1) par->distort=DISTORT;
  if (par->filetype==-2) par->filetype=FILETYPE;
  if (par->ftyplya==-1) par->ftyplya=FTYPLYA;
  if (par->ftypred==-1) par->ftypred=FTYPRED;
  if (par->helio==-1) par->helio=HELIO;
  if (par->linear==-1) par->linear=LINEAR;
  if (par->macmap==-1) par->macmap=MACMAP;
  if (par->nocont==-1) par->nocont=NOCONT;
  if (par->nordclip==-1) par->nordclip=NORDCLIP;
  if (par->nordsig==-1) par->nordsig=NORDSIG;
  if (par->rankspec==-1) par->rankspec=RANKSPEC;
  if (par->raw==-1) par->raw=RAW;
  if (par->replay==-1) par->replay=REPLAY;
  if (par->save==-1) par->save=SAVE;
  if (par->scalmeth==-1) par->scalmeth=SCALMETH;
  if (par->nscalclip==-1) par->nscalclip=NSCALCLIP;
  if (par->thar==-1) par->thar=THAR;
  if (par->vacwl==-1) par->vacwl=VACWL;
  if (par->vshift==-1) par->vshift=VSHIFT;
  if (par->backvers==-1) par->backvers=0;

  if (!par->macmapfile) sprintf(par->macmapfile,"%s",MACMAPFILE);
  if (!par->vshiftfile) sprintf(par->vshiftfile,"%s",VSHIFTFILE);

  return 1;

}
