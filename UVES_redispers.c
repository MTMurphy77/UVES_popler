/****************************************************************************
* Redispers all orders in a spectrum according to the vac-helio wavelength
* scale of the combined spectrum
*
* Note: ranseed is the random number seed which determines the
* "offset" and "slope" parameters of the velocity distortions, when
* requested by the user, in UVES_wpol. It should always come from the
* same spectrum (e.g. the first spectrum read in), but this can
* obviously be changed however the calling routine sees fit. Just be
* sure that the routine which calls UVES_redispers uses the same
* spectrum as those which call UVES_wpol or UVES_revwpol directly.
****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

/* TEMPORARY: The following to be used when creating ASCII spectra of
   individual orders */
/*
#include "file.h"
*/

#define ERR_ARRAY { \
    errormsg("UVES_redispers(): Cannot allocate memory for\n\
\t%s array for redispersion of size %d for\n\
\tspectrum in file\n\t%s",array_name,array_size,filename); }
#define NFERR_ARRAY { \
    nferrormsg("UVES_memspec(): Cannot allocate memory for\n\
\t%s array of size %d\n\tfor order %d for spectrum in file\n\t%s",\
	       array_name,array_size,l+1,filename); \
    return 0; }

int UVES_redispers(spectrum *spec, cspectrum *cspec, long ranseed, params *par) {

  double   wgt=0.0,swgt=0.0,sfrac=0.0,swgtf=0.0,sfracr=0.0;
  double   sfracf=0.0,sfracesq=0.0,sfracesqi=0.0,sfracsqesq=0.0;
  double   swl=0.0,ewl=0.0,rdswl=0.0,rdewl=0.0,lwl=0.0,rwl=0.0;
  double   fr=0.0,sfr=0.0,efr=0.0;
  double   *dat=NULL,*err=NULL,*res=NULL,*frac=NULL; /* Data arrays for averaging */
  double   *cdat=NULL,*cerr=NULL,*cres=NULL,*cfrac=NULL; /* and sig-clipped values */
  double   *tdat=NULL,*terr=NULL,*tfrac=NULL; /* ThAr data arrays for averaging */
  double   *ctdat=NULL,*cterr=NULL,*ctfrac=NULL; /* ThAr sig-clipped values */
  int      array_size=0;
  int      fval=0,lval=0; /* First and last valid pixels in raw order */
  int      ndat=0,ncdat=0,cst=0,maxndat=0;
  int      ntdat=0,nctdat=0,ctst=0;
  int      sidx=0,eidx=0;
  int      i=0,j=0,k=0,l=0;
  char     array_name[NAMELEN]="\0",filename[LNGSTRLEN]="\0";

  /* TEMPORARY: The following to be used when creating ASCII spectra
     of individual orders */
  /*
  char     infile[LNGSTRLEN]="\0";
  char     *cptr;
  FILE     *data_file=NULL;
  */

  /* Check that cspectrum's wavelength scale has been set */
  if (cspec->wl[0]<DRNDTOL)
    errormsg("UVES_redispers(): Combined spectrum's wavelength\n\
\tscale has not been set");

  /** Object spectrum **/
  /* Loop over all orders */
  for (l=0; l<spec->nor; l++) {
    /* Check first to see if order is useful */
    if (spec->or[l].nuse>=MINUSE) {
      /* Find first valid pixel in raw array */
      i=0; while (spec->or[l].st[i]<1) i++; fval=i;
      /* Find last valid pixel in raw array */
      i=spec->or[l].np-1; while (spec->or[l].st[i]<1) i--; lval=i;

      /* TEMPORARY: Create an ASCII file for each echelle order */
      /*
      cptr=strstr(spec->abfile,"fxb_"); cptr+=4; sprintf(infile,"%s",cptr);
      // cptr=strstr(spec->abfile,"HARPN."); cptr+=6; sprintf(infile,"%s",cptr);
      cptr=strstr(infile,".fits"); sprintf(cptr,"_%02d.dat",l+1);
      data_file=faskwopen("File name for echelle order?",infile,4);
      for (i=fval; i<=lval; i++)
	fprintf(data_file,"%-4d  %14.8lf  %14.6lf  %14.6lf\n",i,spec->or[l].vhwl[i],
		spec->or[l].fl[i],spec->or[l].er[i]);
      fclose(data_file);
      */

      /* Find first pixel in combined array to which first valid pixel
	 in raw array will contribute */
      swl=(fval) ? spec->or[l].vhrwl[fval-1] : spec->or[l].fvhlwl;
      if ((spec->or[l].csidx=idxdval(cspec->rwl,cspec->np,swl))==-1)
	errormsg("UVES_redispers(): Cannot find wavelength\n\
\t%9.4lf (left edge of vachelio pixel %d in order %d of file\n\t%s)\n\
\tin combined spectrum wavelength array covering range\n\
\t%9.4lf to %9.4lf",swl,fval+1,l+1,spec->file,cspec->wl[0],
		 cspec->wl[cspec->np-1]);

      /* Find last pixel in combined array to which last pixel in raw array
	 will contribute */
      ewl=spec->or[l].vhrwl[lval];
      if ((spec->or[l].ceidx=idxdval(cspec->rwl,cspec->np,ewl))==-1)
	errormsg("UVES_redispers(): Cannot find wavelength\n\
\t%9.4lf (right edge of vachelio pixel %d in order %d of file\n\t%s)\n\
\tin combined spectrum wavelength array covering range\n\
\t%9.4lf to %9.4lf",ewl,lval+1,l+1,spec->file,cspec->wl[0],
		 cspec->wl[cspec->np-1]);
      if (ewl*(1.0-DRNDTOL)<=cspec->rwl[spec->or[l].ceidx-1]) spec->or[l].ceidx--;
      
      /* Find number of cspec pixels covered by order */
      if ((spec->or[l].nrdp=spec->or[l].ceidx-spec->or[l].csidx+1)<=0) {
	nferrormsg("UVES_redispers(): Found strange number of\n\
\tcombined spectrum pixels (=%d) traversed by valid pixels in\n\
\torder %d of file\n\t%s",spec->or[l].nrdp,l+1,spec->file); return 0;
      }

      /* Allocate memory to redispersed arrays */
      array_size=spec->or[l].nrdp; sprintf(filename,"%s",spec->file);
      sprintf(array_name,"redispersed flux");
      if ((spec->or[l].rdfl=darray(array_size))==NULL) { NFERR_ARRAY; }
      sprintf(array_name,"redispersed error");
      if ((spec->or[l].rder=darray(array_size))==NULL) { NFERR_ARRAY; }
      sprintf(array_name,"redispersed expected fluctuation");
      if ((spec->or[l].rdef=darray(array_size))==NULL) { NFERR_ARRAY; }
      /* "Alternative" arrays are for pre-v0.74 backwards compatibility */
      if (par->version<0.74) {
	sprintf(array_name,"alt. redispersed flux");
	if ((spec->or[l].rdfl_a=darray(array_size))==NULL) { NFERR_ARRAY; }
	sprintf(array_name,"alt. redispersed error");
	if ((spec->or[l].rder_a=darray(array_size))==NULL) { NFERR_ARRAY; }
	sprintf(array_name,"alt. redispersed expected fluctuation");
	if ((spec->or[l].rdef_a=darray(array_size))==NULL) { NFERR_ARRAY; }
      }
      sprintf(array_name,"redispersed continuum");
      if ((spec->or[l].rdco=darray(array_size))==NULL) { NFERR_ARRAY; }
      sprintf(array_name,"old redispersed continuum");
      if ((spec->or[l].ordco=darray(array_size))==NULL) { NFERR_ARRAY; }
      sprintf(array_name,"redispersed median");
      if ((spec->or[l].rdme=darray(array_size))==NULL) { NFERR_ARRAY; }
      /* "Alternative" arrays are for pre-v0.74 backwards compatibility */
      if (par->version<0.74) {
	sprintf(array_name,"alt. redispersed median");
	if ((spec->or[l].rdme_a=darray(array_size))==NULL) { NFERR_ARRAY; }
      }
      sprintf(array_name,"redispersed resolution");
      if ((spec->or[l].rdres=darray(array_size))==NULL) { NFERR_ARRAY; }
      if (par->thar==1) {
	sprintf(array_name,"redispersed ThAr flux");
	if ((spec->or[l].rdth=darray(array_size))==NULL) { NFERR_ARRAY; }
	sprintf(array_name,"redispersed ThAr error");
	if ((spec->or[l].rdter=darray(array_size))==NULL) { NFERR_ARRAY; }
      }
      sprintf(array_name,"redispersed status");
      if ((spec->or[l].rdst=iarray(array_size))==NULL) { NFERR_ARRAY; }
      sprintf(array_name,"old redispersed status");
      if ((spec->or[l].ordst=iarray(array_size))==NULL) { NFERR_ARRAY; }
      if (par->thar==1) {
	sprintf(array_name,"redispersed ThAr status");
	if ((spec->or[l].rdtst=iarray(array_size))==NULL) { NFERR_ARRAY; }
      }

      /* Find maximum size of data arrays needed and allocate memory */
      maxndat=10;
      if (spec->or[l].vhrwl[0]-spec->or[l].fvhlwl>0.0)
	maxndat=(int)((cspec->wl[spec->or[l].csidx+1]-cspec->wl[spec->or[l].csidx])/
		      (spec->or[l].vhrwl[0]-spec->or[l].fvhlwl));
      if (spec->or[l].vhrwl[spec->or[l].np-1]-spec->or[l].vhrwl[spec->or[l].np-2]>0.0)
	maxndat=(int)(MAX(maxndat,
		((cspec->wl[spec->or[l].ceidx]-cspec->wl[spec->or[l].ceidx-1])/
		 (spec->or[l].vhrwl[spec->or[l].np-1]-
		  spec->or[l].vhrwl[spec->or[l].np-2]))))+10;
      maxndat=array_size=(MAX(maxndat,10));
      sprintf(array_name,"dat"); if ((dat=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"err"); if ((err=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"res"); if ((res=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"frac"); if ((frac=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"cdat"); if ((cdat=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"cerr"); if ((cerr=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"cres"); if ((cres=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"cfrac"); if ((cfrac=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"tdat"); if ((tdat=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"terr"); if ((terr=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"tfrac"); if ((tfrac=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"ctdat"); if ((ctdat=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"cterr"); if ((cterr=darray(maxndat))==NULL) { ERR_ARRAY; }
      sprintf(array_name,"ctfrac"); if ((ctfrac=darray(maxndat))==NULL) { ERR_ARRAY; }

      /** Do the redispersion **/
      /* Loop over redispersed pixels */
      for (i=0,j=spec->or[l].csidx; i<spec->or[l].nrdp; i++,j++) {
	/* Find the first and last raw pixels which contribute to this
	   redispersed pixel */
	rdswl=(j) ? cspec->rwl[j-1] : cspec->flwl; rdewl=cspec->rwl[j];
	sidx=(i) ? eidx : fval;
	if (i==spec->or[l].nrdp-1) eidx=lval;
	else {
	  if ((eidx=idxdval(&(spec->or[l].vhrwl[sidx]),spec->or[l].np-sidx,rdewl))==-1)
	    errormsg("UVES_redispers(): Cannot find wavelength %9.4lf\n\
\t(right edge of redispersed pixel %d) in raw vac-helio array\n\
\tranging over %9.4lf to %9.4lf",rdewl,j+1,spec->or[l].vhrwl[sidx],
		   spec->or[l].vhrwl[spec->or[l].np-1]);
	  eidx+=sidx;
	}

	/* Loop through relevant raw pixels and find weighted mean flux */
	ndat=ncdat=ntdat=nctdat=0; cst=ctst=RCLIP; for (k=sidx; k<=eidx; k++) {
	  /* Find fraction of pixel (in pixel space) within starting
	     and ending wavelengths */
	  lwl=(k) ? spec->or[l].vhrwl[k-1] : spec->or[l].fvhlwl;
	  if (lwl>=rdswl) sfr=0.0;
	  else sfr=0.5-(double)k+UVES_revwpol(spec,l,k,rdswl,ranseed,0,par);
	  if (sfr<0.0 || sfr>1.0)
	    errormsg("UVES_redisperse(): Fraction of raw pixel %d\n\
\tfalling bluewards of combined pixel %d in order %d of spectrum %d,\n\t%s,\n\
\tis %lg. Must be >=0.0 and <=1.0",k+1,j+1,l+1,spec->id+1,spec->file,sfr); 
	  rwl=spec->or[l].vhrwl[k];
	  if (rwl<rdewl) efr=0.0;
	  else efr=(double)k+0.5-UVES_revwpol(spec,l,k,rdewl,ranseed,0,par);
	  if (efr<0.0 || efr>1.0)
	    errormsg("UVES_redisperse(): Fraction of raw pixel %d\n\
\tfalling redwards of combined pixel %d in order %d of spectrum %d,\n\t%s,\n\
\tis %lg. Must be >=0.0 and <=1.0",k+1,j+1,l+1,spec->id+1,spec->file,efr); 
	  if ((fr=1.0-sfr-efr)<0.0 || fr>1.0)
	    errormsg("UVES_redisperse(): Fraction of raw pixel %d\n\
\tin order %d of spectrum %d,\n\t%s,\n\tcovering combined pixel %d is %lg.\n\
\tMust be >=0.0 and <=1.0.",k+1,l+1,spec->id+1,spec->file,j+1,fr);
	  fr=(MAX(fr,DRNDTOL));
	  if (spec->or[l].st[k]==1) {
	    dat[ndat]=spec->or[l].fl[k]; err[ndat]=spec->or[l].er[k];
	    if (spec->ftype==FTUVES && spec->fvers==1) res[ndat]=spec->or[l].res[k];
	    frac[ndat++]=fr;
	  } else {
	    if (spec->or[l].st[k]==OCLIP || spec->or[l].st[k]==ACLIP) {
	      if (cst==RCLIP) cst=spec->or[l].st[k];
	      cdat[ncdat]=spec->or[l].fl[k]; cerr[ncdat]=spec->or[l].er[k];
	      if (spec->ftype==FTUVES && spec->fvers==1) res[ncdat]=spec->or[l].res[k];
	      cfrac[ncdat++]=fr;
	    }
	  }
	  if (par->thar==1) {
	    if (spec->or[l].tst[k]==1) {
	      tdat[ntdat]=spec->or[l].th[k]; terr[ntdat]=spec->or[l].ter[k];
	      tfrac[ntdat++]=fr;
	    } else {
	      if (spec->or[l].tst[k]==OCLIP || spec->or[l].st[k]==ACLIP) {
		if (ctst==RCLIP) ctst=spec->or[l].tst[k];
		ctdat[nctdat]=spec->or[l].th[k]; cterr[nctdat]=spec->or[l].ter[k];
		ctfrac[nctdat++]=fr;
	      }
	    }
	  }
	}

	/* Compute statistics for object flux and error */
	sfrac=swgt=swgtf=sfracr=0.0;
	if (!ndat) {
	  if (!ncdat) {
	    spec->or[l].rdfl[i]=0.0;
	    spec->or[l].rder[i]=spec->or[l].rdef[i]=spec->or[l].rdres[i]=-INFIN;
	  } else {
	    for (k=0; k<ncdat; k++) {
	      sfrac+=cfrac[k]; swgt+=(wgt=cfrac[k]/cerr[k]/cerr[k]);
	      swgtf+=wgt*cdat[k];
	      if (spec->ftype==FTUVES && spec->fvers==1) sfracr+=cfrac[k]*cres[k];
	    }
	    spec->or[l].rdfl[i]=swgtf/swgt; spec->or[l].rder[i]=sqrt(1.0/swgt);
	    spec->or[l].rdef[i]=sqrt(sfrac/swgt);
	  }
	  if (spec->ftype==FTUVES && spec->fvers==1)
	    spec->or[l].rdres[i]=sfracr/sfrac;
	  else spec->or[l].rdres[i]=-INFIN;
	  spec->or[l].rdst[i]=cst;
	  if (spec->or[l].rdst[i]!=RCLIP && spec->or[l].rdst[i]!=OCLIP &&
	      spec->or[l].rdst[i]!=ACLIP)
	    errormsg("UVES_redispers(): Strange value of average\n\
\tstatus (=%d) for redispersed pixel %d of order %d.\n\
\tShould be either %d, %d or %d",spec->or[l].rdst[i],i+1,l+1,RCLIP,OCLIP,ACLIP);
	} else {
	  for (k=0; k<ndat; k++) {
	    sfrac+=frac[k]; swgt+=(wgt=frac[k]/err[k]/err[k]); swgtf+=wgt*dat[k];
	    if (spec->ftype==FTUVES && spec->fvers==1) sfracr+=frac[k]*res[k];
	  }
	  spec->or[l].rdfl[i]=swgtf/swgt; spec->or[l].rder[i]=sqrt(1.0/swgt);
	  spec->or[l].rdef[i]=sqrt(sfrac/swgt); spec->or[l].rdst[i]=1;
	  if (spec->ftype==FTUVES && spec->fvers==1) spec->or[l].rdres[i]=sfracr/sfrac;
	  else spec->or[l].rdres[i]=-INFIN;
	}

	/* Compute statistics for alternative object flux and error
	   for pre-0.74 backwards compatibility */
	if (par->version<0.74) {
	  spec->or[l].rdfl_a[i]=spec->or[l].rdfl[i];
	  spec->or[l].rder_a[i]=spec->or[l].rder[i];
	  spec->or[l].rdef_a[i]=spec->or[l].rdef[i];
	  sfrac=sfracf=sfracesq=sfracsqesq=0.0;
	  if (!ndat) {
	    if (!ncdat) {
	      spec->or[l].rdfl[i]=0.0;
	      spec->or[l].rder[i]=spec->or[l].rdef[i]=-INFIN;
	    } else {
	      for (k=0; k<ncdat; k++) {
		sfrac+=cfrac[k]; sfracf+=cfrac[k]*cdat[k];
		sfracesq+=(sfracesqi=cfrac[k]*cerr[k]*cerr[k]);
		sfracsqesq+=cfrac[k]*sfracesqi;
	      }
	      spec->or[l].rdfl[i]=sfracf/sfrac;
	      spec->or[l].rder[i]=sqrt(sfracesq/sfrac);
	      spec->or[l].rdef[i]=sqrt(sfracsqesq)/sfrac;
	    }
	  } else {
	    for (k=0; k<ndat; k++) {
	      sfrac+=frac[k]; sfracf+=frac[k]*dat[k];
	      sfracesq+=(sfracesqi=frac[k]*err[k]*err[k]);
	      sfracsqesq+=frac[k]*sfracesqi;
	    }
	    spec->or[l].rdfl[i]=sfracf/sfrac;
	    spec->or[l].rder[i]=sqrt(sfracesq/sfrac);
	    spec->or[l].rdef[i]=sqrt(sfracsqesq)/sfrac;
	  }
	}

	/* Compute statistics for ThAr flux and error */
	if (par->thar==1) {
	  sfrac=sfracf=sfracesq=0.0;
	  if (!ntdat) {
	    if (!nctdat) {
	      spec->or[l].rdth[i]=0.0; spec->or[l].rdter[i]=-INFIN;
	    } else {
	      for (k=0; k<nctdat; k++) {
		sfrac+=ctfrac[k]; swgt+=(wgt=ctfrac[k]/cterr[k]/cterr[k]);
		swgtf+=wgt*ctdat[k];
	      }
	      spec->or[l].rdth[i]=swgtf/swgt; spec->or[l].rdter[i]=sqrt(1.0/swgt);
	    }
	    spec->or[l].rdtst[i]=ctst;
	    if (spec->or[l].rdtst[i]!=RCLIP && spec->or[l].rdtst[i]!=OCLIP &&
		spec->or[l].rdtst[i]!=ACLIP)
	      errormsg("UVES_redispers(): Strange value of average\n\
\tstatus (=%d) for redispersed ThAr pixel %d of order %d.\n\
\tShould be either %d, %d or %d",spec->or[l].rdtst[i],i+1,l+1,RCLIP,OCLIP,ACLIP);
	  } else {
	    for (k=0; k<ntdat; k++) {
	      sfrac+=tfrac[k]; swgt+=(wgt=tfrac[k]/terr[k]/terr[k]);
	      swgtf+=wgt*tdat[k];
	    }
	    spec->or[l].rdth[i]=swgtf/swgt;
	    spec->or[l].rdter[i]=sqrt(1.0/swgt);
	    spec->or[l].rdtst[i]=1;
	  }
	}
      }

      /* Clean up temporary arrays */
      free(dat); free(err); free(res); free(frac); free(cdat); free(cerr); free(cres);
      free(cfrac); free(tdat); free(terr); free(tfrac); free(ctdat); free(cterr);
      free(ctfrac);

      /* Free space taken by original data arrays (only redispersed values
	 used from here on) */
      if (!UVES_memspec(spec,par,-l-1,0))
	errormsg("UVES_redispers(): Error returned from UVES_memspec() when\n\
\tattempting to free memory from order %d of file\n\t%s",l+1,spec->file);
    }
  }

  /** Sky spectrum **/
  if (spec->skysub) {
    /* Loop over all orders */
    for (l=0; l<spec->ss.nor; l++) {
      /* Find first valid pixel in raw sky order whose right-hand edge
	 wavelength lies redwards of the left-hand edge of first
	 combined pixel */
      i=0;
      while (i<spec->ss.or[l].np && spec->ss.or[l].st[i]<1 &&
	     spec->ss.or[l].vhrwl[i]<cspec->flwl) i++;
      fval=i;
      /* Find last valid pixel in raw sky order whose left-hand edge
	 wavelength lies bluewards of the right-hand edge of last
	 combined pixel */
      i=spec->ss.or[l].np-1;
      while (i>=0 && spec->ss.or[l].st[i]<1 &&
	     spec->ss.or[l].vhrwl[i-1]>cspec->rwl[cspec->np-1]) i--;
      lval=i;
      /* Do not continue with redispersion if all of the sky order
	 lies completely outside the combined spectrum wavelength
	 range */
      if (fval==spec->ss.or[l].np || lval==-1) spec->ss.or[l].nuse=0;
      else spec->ss.or[l].nuse=lval-fval+1;
      if (spec->ss.or[l].nuse) {
	/* Find first pixel in combined array to which first valid pixel
	   in raw sky array will contribute */
	swl=(fval) ? spec->ss.or[l].vhrwl[fval-1] : spec->ss.or[l].fvhlwl;
	if ((spec->ss.or[l].csidx=idxdval(cspec->rwl,cspec->np,swl))==-1)
	  errormsg("UVES_redispers(): Cannot find wavelength\n\
\t%9.4lf (left edge of vachelio sky pixel %d in order %d of file\n\t%s)\n\
\tin combined spectrum wavelength array covering range\n\
\t%9.4lf to %9.4lf",swl,fval+1,l+1,spec->ss.file,cspec->wl[0],
		   cspec->wl[cspec->np-1]);
	/* Find last pixel in combined array to which last valid pixel
	   in raw sky array will contribute */
	ewl=spec->ss.or[l].vhrwl[lval];
	if ((spec->ss.or[l].ceidx=idxdval(cspec->rwl,cspec->np,ewl))==-1)
	  spec->ss.or[l].ceidx=cspec->np-1;
	  /*
	  errormsg("UVES_redispers(): Cannot find wavelength\n\
\t%9.4lf (right edge of vachelio sky pixel %d in order %d of file\n\t%s)\n\
\tin combined spectrum wavelength array covering range\n\
\t%9.4lf to %9.4lf",ewl,lval+1,l+1,spec->ss.file,cspec->wl[0],
		   cspec->wl[cspec->np-1]);
	  */
	if (ewl*(1.0-DRNDTOL)<=cspec->rwl[spec->ss.or[l].ceidx-1])
	  spec->ss.or[l].ceidx--;
	/* Find number of cspec pixels covered by sky order */
	if ((spec->ss.or[l].nrdp=spec->ss.or[l].ceidx-spec->ss.or[l].csidx+1)<=0) {
	  nferrormsg("UVES_redispers(): Found strange number of\n\
\tcombined spectrum pixels (=%d) traversed by valid sky pixels in\n\
\torder %d of file\n\t%s",spec->ss.or[l].nrdp,l+1,spec->ss.file); return 0;
	}
	/* Allocate memory to redispersed sky arrays */
	array_size=spec->ss.or[l].nrdp; sprintf(filename,"%s",spec->ss.file);
	sprintf(array_name,"redispersed sky flux");
	if ((spec->ss.or[l].rdfl=darray(array_size))==NULL) { NFERR_ARRAY; }
	sprintf(array_name,"redispersed sky error");
	if ((spec->ss.or[l].rder=darray(array_size))==NULL) { NFERR_ARRAY; }
	sprintf(array_name,"redispersed sky expected fluctuation");
	if ((spec->ss.or[l].rdef=darray(array_size))==NULL) { NFERR_ARRAY; }
	sprintf(array_name,"redispersed sky median");
	if ((spec->ss.or[l].rdme=darray(array_size))==NULL) { NFERR_ARRAY; }
	sprintf(array_name,"redispersed sky status");
	if ((spec->ss.or[l].rdst=iarray(array_size))==NULL) { NFERR_ARRAY; }
	/* Find maximum size of data arrays needed and allocate memory */
	maxndat=10;
	if (spec->ss.or[l].vhrwl[0]-spec->ss.or[l].fvhlwl>0.0)
	  maxndat=(int)((cspec->wl[spec->ss.or[l].csidx+1]-cspec->wl[spec->ss.or[l].csidx])/
			(spec->ss.or[l].vhrwl[0]-spec->ss.or[l].fvhlwl));
	if (spec->ss.or[l].vhrwl[spec->ss.or[l].np-1]-
	    spec->ss.or[l].vhrwl[spec->ss.or[l].np-2]>0.0)
	  maxndat=(int)(MAX(maxndat,
			    ((cspec->wl[spec->ss.or[l].ceidx]-
			      cspec->wl[spec->ss.or[l].ceidx-1])/
			     (spec->ss.or[l].vhrwl[spec->ss.or[l].np-1]-
			      spec->ss.or[l].vhrwl[spec->ss.or[l].np-2]))))+10;
	maxndat=array_size=(MAX(maxndat,10));
	sprintf(array_name,"sky dat"); if ((dat=darray(maxndat))==NULL) { ERR_ARRAY; }
	sprintf(array_name,"sky err"); if ((err=darray(maxndat))==NULL) { ERR_ARRAY; }
	sprintf(array_name,"sky frac"); if ((frac=darray(maxndat))==NULL) { ERR_ARRAY; }
	sprintf(array_name,"sky cdat"); if ((cdat=darray(maxndat))==NULL) { ERR_ARRAY; }
	sprintf(array_name,"sky cerr"); if ((cerr=darray(maxndat))==NULL) { ERR_ARRAY; }
	sprintf(array_name,"sky cfrac"); if ((cfrac=darray(maxndat))==NULL) { ERR_ARRAY; }
	/** Do the redispersion **/
	/* Loop over redispersed pixels */
	for (i=0,j=spec->ss.or[l].csidx; i<spec->ss.or[l].nrdp; i++,j++) {
	  /* Find the first and last raw pixels which contribute to this
	     redispersed pixel */
	  rdswl=(j) ? cspec->rwl[j-1] : cspec->flwl; rdewl=cspec->rwl[j];
	  sidx=(i) ? eidx : fval;
	  if (i==spec->ss.or[l].nrdp-1) eidx=lval;
	  else {
	    if ((eidx=idxdval(&(spec->ss.or[l].vhrwl[sidx]),spec->ss.or[l].np-sidx,rdewl))==-1)
	      errormsg("UVES_redispers(): Cannot find wavelength %9.4lf\n\
\t(right edge of redispersed sky pixel %d) in raw vac-helio array\n\
\tranging over %9.4lf to %9.4lf",rdewl,j+1,spec->ss.or[l].vhrwl[sidx],
		       spec->ss.or[l].vhrwl[spec->ss.or[l].np-1]);
	    eidx+=sidx;
	  }
	  /* Loop through relevant raw pixels and find weighted mean flux */
	  ndat=ncdat=0; cst=RCLIP; for (k=sidx; k<=eidx; k++) {
	    /* Find fraction of pixel (in pixel space) within starting
	       and ending wavelengths */
	    lwl=(k) ? spec->ss.or[l].vhrwl[k-1] : spec->ss.or[l].fvhlwl;
	    if (lwl>=rdswl) sfr=0.0;
	    else sfr=0.5-(double)k+UVES_revwpol(spec,l,k,rdswl,ranseed,1,par);
	    if (sfr<0.0 || sfr>1.0)
	      errormsg("UVES_redisperse(): Fraction of raw sky pixel %d\n\
\tfalling bluewards of combined pixel %d in order %d of spectrum %d,\n\t%s,\n\
\tis %lg. Must be >=0.0 and <=1.0",k+1,j+1,l+1,spec->ss.or[l].id+1,spec->ss.file,sfr);
	    rwl=spec->ss.or[l].vhrwl[k];
	    if (rwl<rdewl) efr=0.0;
	    else efr=(double)k+0.5-UVES_revwpol(spec,l,k,rdewl,ranseed,1,par);
	    if (efr<0.0 || efr>1.0)
	      errormsg("UVES_redisperse(): Fraction of raw sky pixel %d\n\
\tfalling redwards of combined pixel %d in order %d of spectrum %d,\n\t%s,\n\
\tis %lg. Must be >=0.0 and <=1.0",k+1,j+1,l+1,spec->ss.or[l].id+1,spec->ss.file,efr); 
	    if ((fr=1.0-sfr-efr)<0.0 || fr>1.0)
	      errormsg("UVES_redisperse(): Fraction of raw sky pixel %d\n\
\tin order %d of spectrum %d,\n\t%s,\n\tcovering combined pixel %d is %lg.\n \
\tMust be >=0.0 and <=1.0.",k+1,l+1,spec->ss.or[l].id+1,spec->ss.file,j+1,fr);
	    fr=(MAX(fr,DRNDTOL));
	    if (spec->ss.or[l].st[k]==1) {
	      dat[ndat]=spec->ss.or[l].fl[k]; err[ndat]=spec->ss.or[l].er[k];
	      frac[ndat++]=fr;
	    } else {
	      if (spec->ss.or[l].st[k]==OCLIP || spec->ss.or[l].st[k]==ACLIP) {
		if (cst==RCLIP) cst=spec->ss.or[l].st[k];
		cdat[ncdat]=spec->ss.or[l].fl[k]; cerr[ncdat]=spec->ss.or[l].er[k];
		cfrac[ncdat++]=fr;
	      }
	    }
	  }
	  /* Compute statistics for object flux and error */
	  sfrac=swgt=swgtf=0.0;
	  if (!ndat) {
	    if (!ncdat) {
	      spec->ss.or[l].rdfl[i]=0.0;
	      spec->ss.or[l].rder[i]=spec->ss.or[l].rdef[i]=-INFIN;
	    } else {
	      for (k=0; k<ncdat; k++) {
		sfrac+=cfrac[k]; swgt+=(wgt=cfrac[k]/cerr[k]/cerr[k]);
		swgtf+=wgt*cdat[k];
	      }
	      spec->ss.or[l].rdfl[i]=swgtf/swgt;
	      spec->ss.or[l].rder[i]=sqrt(1.0/swgt);
	      spec->ss.or[l].rdef[i]=sqrt(sfrac/swgt);
	    }
	    spec->ss.or[l].rdst[i]=cst;
	    if (spec->ss.or[l].rdst[i]!=RCLIP && spec->ss.or[l].rdst[i]!=OCLIP &&
		spec->ss.or[l].rdst[i]!=ACLIP)
	      errormsg("UVES_redispers(): Strange value of average\n\
\tstatus (=%d) for redispersed sky pixel %d of order %d.\n\
\tShould be either %d, %d or %d",spec->ss.or[l].rdst[i],i+1,l+1,RCLIP,OCLIP,ACLIP);
	  } else {
	    for (k=0; k<ndat; k++) {
	      sfrac+=frac[k]; swgt+=(wgt=frac[k]/err[k]/err[k]); swgtf+=wgt*dat[k];
	    }
	    spec->ss.or[l].rdfl[i]=swgtf/swgt;
	    spec->ss.or[l].rder[i]=sqrt(1.0/swgt);
	    spec->ss.or[l].rdef[i]=sqrt(sfrac/swgt); spec->ss.or[l].rdst[i]=1;
	  }
	}
	/* Clean up temporary arrays */
	free(dat); free(err); free(frac); free(cdat); free(cerr); free(cfrac);
      }
      /* Free space taken by original data arrays (only redispersed values
	 used from here on) */
      if (!UVES_memspec(spec,par,-l-1,1))
	errormsg("UVES_redispers(): Error returned from UVES_memspec() when\n\
\tattempting to free memory from sky order %d of file\n\t%s",l+1,spec->ss.file);
    }
  }

  return 1;

}
