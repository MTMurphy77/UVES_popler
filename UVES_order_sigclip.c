/****************************************************************************
* Sigma-clip positive points (and neighbours) from an order
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "stats.h"
#include "memory.h"
#include "error.h"

int UVES_order_sigclip(echorder *ord, params *par) {

  double   max=0.0;
  double   *data=NULL;    /* Array on which statistics are determined */
  int      fval=0,lval=0,nval=0,noson2=0;
  int      i=0,j=0,k=0,l=0;
  int      *clip=NULL,*sts=NULL;
  statset  stat;

  /* Check first to see if order is useful */
  if (ord->nuse<MINUSE) return 1;

  /* Allocate memory to data array */
  if ((data=darray(par->nordsig-1))==NULL)
    errormsg("UVES_order_sigclip(): Could not allocate memory\n\
\tfor data array of size %d",par->nordsig-1);
  /* Allocate memory to data array */
  if ((clip=iarray(par->nordsig-1))==NULL)
    errormsg("UVES_order_sigclip(): Could not allocate memory\n\
\tfor clip array of size %d",par->nordsig-1);
  /* Allocate memory to temporary status array and copy order's status array */
  if ((sts=iarray(ord->np))==NULL)
    errormsg("UVES_order_sigclip(): Could not allocate memory\n\
\tfor sts array of size %d",ord->np);
  for (i=0; i<ord->np; i++) sts[i]=ord->st[i];

  /* Find first and last pixels for which statistics can be done */
  noson2=par->nordsig/2;
  for (i=0,fval=0,nval=0; i<ord->np-noson2; i++) {
    if (nval==noson2 && sts[i]>0) { fval=i; break; }
    if (sts[i]>0) nval++;
  }
  if (!fval) errormsg("UVES_order_sigclip(): Cannot find valid starting\n\
\tpixel for sigma-clip of order. Try decreasing nordsig or increasing\n\
\tMINUSE in UVES_popler.h");
  for (i=ord->np-1,lval=0,nval=0; i>noson2; i--) {
    if (nval==noson2 && sts[i]>0) { lval=i; break; }
    if (sts[i]>0) nval++;
  }
  if (!lval) errormsg("UVES_order_sigclip(): Cannot find valid ending\n\
\tpixel for sigma-clip of order. Try decreasing nordsig or increasing\n\
\tMINUSE in UVES_popler.h");
  /* Check to make sure that there are enough valid pixels either side
     of the first and last pixels for which statistics look possible */
  if (fval>=lval) errormsg("UVES_order_sigclip(): Cannot find enough valid\n\
\tpixels in order. Try decreasing nordsig or increasing MINUSE in UVES_popler.h");

  /** Run though order and find positive-deviant points **/
  for (i=fval; i<=lval; i++) {
    /* For each central pixel, find array of nordsig valid pixels */
    for (j=i-1,nval=0; j>=0; j--) {
      if (sts[j]>0) { data[nval]=ord->fl[j]; clip[nval++]=1; }
      if (nval==noson2) break;
    }
    if (j<0) errormsg("UVES_order_sigclip(): Cannot find correct number of\n\
\tvalid pixels to left of first pixel for which clipping is possible.\n\
\tCounting error in code!");
    for (j=i+1; j<ord->np; j++) {
      if (sts[j]>0) { data[nval]=ord->fl[j]; clip[nval++]=1; }
      if (nval==par->nordsig-1) break;
    }
    if (j==ord->np)
      errormsg("UVES_order_sigclip(): Cannot find correct number of\n\
\tvalid pixels to right of first pixel for which clipping is possible.\n\
\tCounting error in code!");
    /* Exclude the nordclip most positive points */
    for (j=0; j<par->nordclip; j++) {
      k=0; while (!clip[k]) k++; max=-INFIN;
      for (l=k; l<par->nordsig-1; l++)
	if (clip[l] && data[l]>max) { max=data[l]; k=l; }
      clip[k]=0;
    }
    /* Take statistics from remaining pixels */
    if (!stats(data,NULL,NULL,NULL,clip,par->nordsig-1,0,&stat)) {
      nferrormsg("UVES_order_sigclip(): Error returned from stats()\n\
\twhen calculating mean and rms of clipped array"); return 0;
    }
    /* Do the sigma-clip, marking the fact using OCLOP in the status array */
    if (stat.mean>stat.rms*par->ordsigzero &&
	ord->fl[i]>stat.mean+stat.rms*par->ordsig) {
      ord->st[i]=OCLIP;
      /* And sigma-clip the immediate neighbours too */
      if (ord->fl[i-1]>stat.mean+stat.rms*par->ordsignbr) ord->st[i-1]=OCLIP;
      if (ord->fl[i+1]>stat.mean+stat.rms*par->ordsignbr) ord->st[i+1]=OCLIP;
    }
  }

  /* Clean up */
  free(data); free(clip); free(sts);
  
  return 1;

}
