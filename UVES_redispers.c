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
#include <math.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

/* TEMPORARY: The following to be used when creating ASCII spectra of
   individual orders */
/*
#include "file.h"
*/

int UVES_redispers(spectrum *spec, cspectrum *cspec, long ranseed, params *par) {

  double   sfrac=0.0,sfracf=0.0,sfracr=0.0,sfracesq=0.0,sfracesqi=0.0,sfracsqesq=0.0;
  double   swl=0.0,ewl=0.0,rdswl=0.0,rdewl=0.0,lwl=0.0,rwl=0.0;
  double   fr=0.0,sfr=0.0,efr=0.0;
  double   *dat=NULL,*err=NULL,*res=NULL,*frac=NULL; /* Data arrays for averaging */
  double   *cdat=NULL,*cerr=NULL,*cres=NULL,*cfrac=NULL; /* and sig-clipped values */
  double   *tdat=NULL,*terr=NULL,*tfrac=NULL; /* ThAr data arrays for averaging */
  double   *ctdat=NULL,*cterr=NULL,*ctfrac=NULL; /* ThAr sig-clipped values */
  int      fval=0,lval=0; /* First and last valid pixels in raw order */
  int      ndat=0,ncdat=0,cst=0,maxndat=0;
  int      ntdat=0,nctdat=0,ctst=0;
  int      sidx=0,eidx=0;
  int      i=0,j=0,k=0,l=0;

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
      if ((spec->or[l].rdfl=darray(spec->or[l].nrdp))==NULL) {
	nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed flux array for order %d of file\n\t%s",l+1,spec->file);
	return 0;
      }
      if ((spec->or[l].rder=darray(spec->or[l].nrdp))==NULL) {
	nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed error array for order %d of file\n\t%s",l+1,spec->file);
	return 0;
      }
      if ((spec->or[l].rdef=darray(spec->or[l].nrdp))==NULL) {
	nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed expected fluctuation array for order %d of file\n\t%s",l+1,
		   spec->file); return 0;
      }
      if ((spec->or[l].rdco=darray(spec->or[l].nrdp))==NULL) {
	nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed continuum array for order %d of file\n\t%s",l+1,spec->file);
	return 0;
      }
      if ((spec->or[l].ordco=darray(spec->or[l].nrdp))==NULL) {
	nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor old redispersed continuum array for order %d of file\n\t%s",l+1,spec->file);
	return 0;
      }
      if ((spec->or[l].rdme=darray(spec->or[l].nrdp))==NULL) {
	nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed median error array for order %d of file\n\t%s",l+1,spec->file);
	return 0;
      }
      if ((spec->or[l].rdres=darray(spec->or[l].nrdp))==NULL) {
	nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed resolution error array for order %d of file\n\t%s",l+1,spec->file);
	return 0;
      }
      if (par->thar==1) {
	if ((spec->or[l].rdth=darray(spec->or[l].nrdp))==NULL) {
	  nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed ThAr flux array for order %d of file\n\t%s",l+1,spec->file);
	  return 0;
	}
	if ((spec->or[l].rdter=darray(spec->or[l].nrdp))==NULL) {
	  nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed ThAr flux array for order %d of file\n\t%s",l+1,spec->file);
	  return 0;
	}
      }
      if ((spec->or[l].rdst=iarray(spec->or[l].nrdp))==NULL) {
	nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed status array for order %d of file\n\t%s",l+1,spec->file);
	return 0;
      }
      if ((spec->or[l].ordst=iarray(spec->or[l].nrdp))==NULL) {
	nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor old redispersed status array for order %d of file\n\t%s",l+1,spec->file);
	return 0;
      }
      if (par->thar==1) {
	if ((spec->or[l].rdtst=iarray(spec->or[l].nrdp))==NULL) {
	  nferrormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor redispersed ThAr status array for order %d of file\n\t%s",l+1,spec->file);
	  return 0;
	}
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
      maxndat=(MAX(maxndat,10));
      if ((dat=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor dat array for redispersion of size %d",maxndat);
      if ((err=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor err array for redispersion of size %d",maxndat);
      if ((res=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor res array for redispersion of size %d",maxndat);
      if ((frac=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor frac array for redispersion of size %d",maxndat);
      if ((cdat=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor cdat array for redispersion of size %d",maxndat);
      if ((cerr=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor cerr array for redispersion of size %d",maxndat);
      if ((cres=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor cres array for redispersion of size %d",maxndat);
      if ((cfrac=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor cfrac array for redispersion of size %d",maxndat);
      if ((tdat=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor tdat array for redispersion of size %d",maxndat);
      if ((terr=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor terr array for redispersion of size %d",maxndat);
      if ((tfrac=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor tfrac array for redispersion of size %d",maxndat);
      if ((ctdat=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor ctdat array for redispersion of size %d",maxndat);
      if ((cterr=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor cterr array for redispersion of size %d",maxndat);
      if ((ctfrac=darray(maxndat))==NULL)
	errormsg("UVES_redispers(): Cannot allocate memory\n\
\tfor ctfrac array for redispersion of size %d",maxndat);

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
\t(right edge of redispersed pixel %d) in raw vac-helio array\n	       \
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
	  else sfr=0.5-(double)k+UVES_revwpol(spec,l,k,rdswl,ranseed,par);
	  if (sfr<0.0 || sfr>1.0)
	    errormsg("UVES_redisperse(): Fraction of raw pixel %d\n\
\tfalling bluewards of combined pixel %d in order %d of spectrum %d,\n\t%s,\n\
\tis %lg. Must be >=0.0 and <=1.0",k+1,j+1,l+1,spec->id+1,spec->file,sfr); 
	  rwl=spec->or[l].vhrwl[k];
	  if (rwl<rdewl) efr=0.0;
	  else efr=(double)k+0.5-UVES_revwpol(spec,l,k,rdewl,ranseed,par);
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
	sfrac=sfracf=sfracesq=sfracsqesq=sfracr=0.0;
	if (!ndat) {
	  if (!ncdat) {
	    spec->or[l].rdfl[i]=0.0;
	    spec->or[l].rder[i]=spec->or[l].rdef[i]=spec->or[l].rdres[i]=-INFIN;
	  } else {
	    for (k=0; k<ncdat; k++) {
	      sfrac+=cfrac[k]; sfracf+=cfrac[k]*cdat[k];
	      sfracesq+=(sfracesqi=cfrac[k]*cerr[k]*cerr[k]);
	      sfracsqesq+=cfrac[k]*sfracesqi;
	      if (spec->ftype==FTUVES && spec->fvers==1) sfracr+=cfrac[k]*cres[k];
	    }
	    spec->or[l].rdfl[i]=sfracf/sfrac;
	    spec->or[l].rder[i]=sqrt(sfracesq/sfrac);
	    spec->or[l].rdef[i]=sqrt(sfracsqesq)/sfrac;
	    if (spec->ftype==FTUVES && spec->fvers==1)
	      spec->or[l].rdres[i]=sfracr/sfrac;
	    else spec->or[l].rdres[i]=-INFIN;
	  }
	  spec->or[l].rdst[i]=cst;
	  if (spec->or[l].rdst[i]!=RCLIP && spec->or[l].rdst[i]!=OCLIP &&
	      spec->or[l].rdst[i]!=ACLIP)
	    errormsg("UVES_redispers(): Strange value of average\n\
\tstatus (=%d) for redispersed pixel %d of order %d.\n\
\tShould be either %d, %d or %d",spec->or[l].rdst[i],i+1,l+1,RCLIP,OCLIP,ACLIP);
	} else {
	  for (k=0; k<ndat; k++) {
	    sfrac+=frac[k]; sfracf+=frac[k]*dat[k];
	    sfracesq+=(sfracesqi=frac[k]*err[k]*err[k]);
	    sfracsqesq+=frac[k]*sfracesqi;
	    if (spec->ftype==FTUVES && spec->fvers==1) sfracr+=frac[k]*res[k];
	  }
	  spec->or[l].rdfl[i]=sfracf/sfrac;
	  spec->or[l].rder[i]=sqrt(sfracesq/sfrac);
	  spec->or[l].rdef[i]=sqrt(sfracsqesq)/sfrac; spec->or[l].rdst[i]=1;
	  if (spec->ftype==FTUVES && spec->fvers==1) spec->or[l].rdres[i]=sfracr/sfrac;
	  else spec->or[l].rdres[i]=-INFIN;
	}

	/* Compute statistics for ThAr flux and error */
	if (par->thar==1) {
	  sfrac=sfracf=sfracesq=0.0;
	  if (!ntdat) {
	    if (!nctdat) {
	      spec->or[l].rdth[i]=0.0; spec->or[l].rdter[i]=-INFIN;
	    } else {
	      for (k=0; k<nctdat; k++) {
		sfrac+=ctfrac[k]; sfracf+=ctfrac[k]*ctdat[k];
		sfracesq+=ctfrac[k]*cterr[k]*cterr[k];
	      }
	      spec->or[l].rdth[i]=sfracf/sfrac;
	      spec->or[l].rdter[i]=sqrt(sfracesq/sfrac);
	    }
	    spec->or[l].rdtst[i]=ctst;
	    if (spec->or[l].rdtst[i]!=RCLIP && spec->or[l].rdtst[i]!=OCLIP &&
		spec->or[l].rdtst[i]!=ACLIP)
	      errormsg("UVES_redispers(): Strange value of average\n\
\tstatus (=%d) for redispersed ThAr pixel %d of order %d.\n		\
\tShould be either %d, %d or %d",spec->or[l].rdtst[i],i+1,l+1,RCLIP,OCLIP,ACLIP);
	  } else {
	    for (k=0; k<ntdat; k++) {
	      sfrac+=tfrac[k]; sfracf+=tfrac[k]*tdat[k];
	      sfracesq+=tfrac[k]*terr[k]*terr[k];
	    }
	    spec->or[l].rdth[i]=sfracf/sfrac;
	    spec->or[l].rdter[i]=sqrt(sfracesq/sfrac); spec->or[l].rdtst[i]=1;
	  }
	}
      }

      /* Clean up temporary arrays */
      free(dat); free(err); free(res); free(frac); free(cdat); free(cerr); free(cres);
      free(cfrac); free(tdat); free(terr); free(tfrac); free(ctdat); free(cterr);
      free(ctfrac);

      /* Free space taken by original data arrays (only redispersed values
	 used from here on) */
      if (!UVES_memspec(spec,par,-l-1))
	errormsg("UVES_redispers(): Error returned from UVES_memspec() when\n\
\tattempting to free memory from order %d of file\n\t%s",l+1,spec->file);
      /*
      free(spec->or[l].vhwl); free(spec->or[l].vhrwl); free(spec->or[l].fl);
      free(spec->or[l].er); free(spec->or[l].res); free(spec->or[l].st);
      if (par->thar==1) {
	free(spec->or[l].th); free(spec->or[l].ter); free(spec->or[l].tst);
      }
      */
    }
  }
      
  return 1;

}
