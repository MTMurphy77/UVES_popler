/****************************************************************************
* Merge the orders for an individual ThAr exposure in a very simply
* step-function way.
****************************************************************************/

#include "UVES_popler.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int UVES_merge_thar(spectrum *spec, cspectrum *cspec, params *par) {

  int      csidx=0,ceidx=0,sidx=0,eidx=0,osidx=0,oeidx=0;
  int      sp=0,ep=0;
  int      i=0,j=0,k=0,l=0;

  /* Initialize starting and ending pixel indices to the max and min
     values from the combined spectrum */
  csidx=cspec->np-1; ceidx=0;

  /** First determine number of redispersed pixels in merged
      spectrum **/
  /* Loop over all orders */
  for (l=0; l<spec->nor; l++) {
    /* Check first to see if order is useful */
    if (spec->or[l].nuse>=MINUSE) {
      csidx=(MIN(spec->or[l].csidx,csidx)); ceidx=(MAX(spec->or[l].ceidx,ceidx));
    }
  }
  spec->th.np=ceidx-csidx+1; spec->th.csidx=csidx; spec->th.ceidx=ceidx;

  /* Set the wavelength of the first and last pixel in the merged
     spectrum, making sure to adjust it for the heliocentric
     correction (which has already been applied by which should now be
     taken off again */
  spec->th.fwl=cspec->wl[csidx]/(1.0+spec->vhel/C_C_K);
  spec->th.lwl=cspec->wl[ceidx]/(1.0+spec->vhel/C_C_K);

  /* Allocate memory to merged ThAr arrays */
  if ((spec->th.fl=darray(spec->th.np))==NULL)
    errormsg("UVES_merge_thar(): Cannot allocate memory\n\
\tto merged ThAr flux array of size %d for file\n\t%s",spec->th.np,spec->file);
  if ((spec->th.er=darray(spec->th.np))==NULL)
    errormsg("UVES_merge_thar(): Cannot allocate memory\n\
\tto merged ThAr error array of size %d for file\n\t%s",spec->th.np,spec->file);
  if ((spec->th.st=iarray(spec->th.np))==NULL)
    errormsg("UVES_merge_thar(): Cannot allocate memory\n\
\tto merged ThAr status array of size %d for file\n\t%s",spec->th.np,spec->file);

  /* Initialize the merged status array */
  for (i=0; i<spec->th.np; i++) spec->th.st[i]=NCLIP;

  /** Form the merged ThAr spectrum. The algorithm is simply to take
      each order in sequence, find it's overlap with the current merged
      spectrum and write the new order to the merged spectrum from
      half-way along the overlap region **/
  /* Find the first useful order */
  l=0; while (l<spec->nor) { if (spec->or[l].nuse>=MINUSE) break; l++; }
  if (l==spec->nor) errormsg("UVES_merge_thar(): There appears to be\n\
\tno useful orders in file\n\t%s",spec->file);
  /* Initialize starting and ending pixel indices */
  sidx=spec->or[l].csidx; eidx=spec->or[l].ceidx;
  /* Write this first spectrum to merged arrays */
  for (j=0,k=sidx-csidx; j<spec->or[l].nrdp; j++,k++) {
    spec->th.fl[k]=spec->or[l].rdth[j]; spec->th.er[k]=spec->or[l].rdter[j];
    spec->th.st[k]=spec->or[l].rdtst[j];
  }
  osidx=sidx; oeidx=eidx;
  /* Loop through other useful orders and merge them into final spectrum */
  for (i=l; i<spec->nor; i++) {
    /* Check first to see if order is useful */
    if (spec->or[i].nuse>=MINUSE) {
      /* Determine overlap with existing merged spectrum */
      sidx=(MIN(osidx,spec->or[i].csidx)); eidx=(MAX(oeidx,spec->or[i].ceidx));
      if (spec->or[i].csidx<osidx && spec->or[i].ceidx>osidx) {
	ep=(osidx+spec->or[i].ceidx)/2-spec->or[i].csidx;
	for (j=0,k=sidx-csidx; j<ep; j++,k++) {
	  spec->th.fl[k]=spec->or[i].rdth[j]; spec->th.er[k]=spec->or[i].rdter[j];
	  spec->th.st[k]=spec->or[i].rdtst[j];
	}
	spec->th.fl[sidx-csidx+ep]=0.0; spec->th.er[sidx-csidx+ep]=-INFIN;
	spec->th.st[sidx-csidx+ep]=SCLIP;
      } else if (spec->or[i].ceidx>oeidx && spec->or[i].csidx<oeidx) {
	sp=(oeidx+spec->or[i].csidx)/2-spec->or[i].csidx;
	for (j=sp+1,k=spec->or[i].csidx-csidx+sp+1; j<spec->or[i].nrdp; j++,k++) {
	  spec->th.fl[k]=spec->or[i].rdth[j]; spec->th.er[k]=spec->or[i].rdter[j];
	  spec->th.st[k]=spec->or[i].rdtst[j];
	}
	spec->th.fl[spec->or[i].csidx-csidx+sp]=0.0;
	spec->th.er[spec->or[i].csidx-csidx+sp]=-INFIN;
	spec->th.st[spec->or[i].csidx-csidx+sp]=SCLIP;
      } else if (spec->or[i].csidx>oeidx || spec->or[i].ceidx<osidx) {
	for (j=0,k=spec->or[i].csidx-csidx; j<spec->or[l].nrdp; j++,k++) {
	  spec->th.fl[k]=spec->or[l].rdth[j]; spec->th.er[k]=spec->or[l].rdter[j];
	  spec->th.st[k]=spec->or[l].rdtst[j];
	}
      }
      osidx=sidx; oeidx=eidx;
    }
  }

  return 1;

}
