/****************************************************************************
* Read in a UVES_popler log file
****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "UVES_popler.h"
#include "file.h"
#include "error.h"

#define READ_DATA_FILE { \
  idx++; if ((cptr=fgets(buffer,LNGSTRLEN,data_file))==NULL) { fclose(data_file); \
    errormsg("UVES_rUPLfile(): Error reading line %d\n\tin file %s",idx,infile); } \
  }
#define ERR_FMT { \
    errormsg("UVES_rUPLfile(): Incorrect format in line %d\n\tof file %s",idx,infile); \
  }

int UVES_rUPLfile(char *infile, spectrum **spec, int *nspec, action **act,
		  int *nact, cspectrum *cspec, params *par) {

  double ddum=0.0;
  int    idx=0;
  int    idum=0;
  int    i=0,j=0,k=0;
  char   buffer[LNGSTRLEN]="\0",newname[LNGSTRLEN]="\0";
  char   stro[LNGSTRLEN]="\0",strn[LNGSTRLEN]="\0";
  char   *cptr,*cp1,*cp2,*cp3,*cp4;
  macmap mmap;
  FILE   *data_file=NULL;

  /** Open input file, see if it's a UPL file or a list of FITS files
      or a single FITS file **/
  if ((data_file=faskropen("Input UPL file name?",infile,5))==NULL)
    errormsg("UVES_rUPLfile(): Can not open file %s",infile);

  /* Read in first line from UPL file */
  READ_DATA_FILE;
  /* Get version number */
  if (sscanf(&(buffer[22]),"%lf",&(par->version))!=1)
    errormsg("UVES_rUPLfile(): Cannot read UPL version number in\n\tfile %s",
	     infile);
  /* Bug in auto-scaling feature was fixed in version 0.40 in
     UVES_rescale_region.c so UPL files before this should be rebadged
     with the same UPL version number for backwards-compatibility. */
  if (par->version<0.40) par->backvers=1;
  /* Exit if the user is using a UPL file which is from a newer
     version of UVES_popler that the one they are currently
     running */
  if (par->version>VERSION)
    errormsg("UVES_rUPLfile(): The UPL file %s is from a newer version\n\
\t(%.2lf) of UVES_popler than the one currently running (%.2lf).\n\
\tPlease download and install latest version of UVES_popler from\n\t%s",infile,
	     par->version,VERSION,WWW);

  /* Get reduction parameters */
  READ_DATA_FILE; READ_DATA_FILE; READ_DATA_FILE;
  if (par->combmeth<0 && sscanf(buffer," combmeth = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->combmeth<0) par->combmeth=idum;
  READ_DATA_FILE;
  if (par->disp<0.0 && sscanf(buffer," disp = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->disp<0.0) par->disp=ddum;
  if (par->version>=0.01) {
    READ_DATA_FILE;
    if (par->filetype==-2 && sscanf(buffer," filetype = %d",&idum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->filetype==-2) {
      par->filetype=idum;
      if (par->filetype<FTMIX || (par->filetype>FTIESI && par->filetype<FTCOMB) ||
	  par->filetype>FTCOMB)
	errormsg("UVES_rUPLfile(): 'filetype' must be >= %d and\n\
<= %d or = %d in UPL file %s",FTMIX,FTIESI,FTCOMB,infile);
    }
  } else par->filetype=FILETYPE;
  if (par->version>=0.31) {
    READ_DATA_FILE;
    if (par->helio<0 && sscanf(buffer," helio = %d",&idum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->helio<0) par->helio=idum; 
    READ_DATA_FILE;
    if (par->vacwl<0 && sscanf(buffer," vacwl = %d",&idum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->vacwl<0) par->vacwl=idum; 
  } else { par->helio=HELIO; par->vacwl=VACWL; }
  READ_DATA_FILE;
  if (par->zem<0.0 && sscanf(buffer," zem = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->zem<0.0) par->zem=ddum;
  READ_DATA_FILE;
  if (par->linear<0 && sscanf(buffer," linear = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->linear<0) par->linear=idum; 
  READ_DATA_FILE;
 if (par->thar<0 && sscanf(buffer," thar = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->thar<0) par->thar=idum;
  READ_DATA_FILE;
  READ_DATA_FILE;
  if (par->nordclip<0 && sscanf(buffer," nordclip = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->nordclip<0) par->nordclip=idum;
  READ_DATA_FILE;
  if (par->nordsig<0 && sscanf(buffer," nordsig = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->nordsig<0) par->nordsig=idum;
  READ_DATA_FILE;
  if (par->ordsig<0.0 && sscanf(buffer," ordsig = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->ordsig<0.0) par->ordsig=ddum;
  READ_DATA_FILE;
  if (par->ordsignbr<0.0 && sscanf(buffer," ordsignbr = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->ordsignbr<0.0) par->ordsignbr=ddum;
  READ_DATA_FILE;
  if (par->ordsigzero<0.0 && sscanf(buffer," ordsigzero = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->ordsigzero<0.0) par->ordsigzero=ddum;
  if (par->version>0.22) {
    READ_DATA_FILE;
    if (par->ordmedfrac<0.0 && sscanf(buffer," ordmedfrac = %lf",&ddum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->ordmedfrac<0.0) par->ordmedfrac=ddum;
    READ_DATA_FILE;
    if (par->ordmedrej<0.0 && sscanf(buffer," ordmedrej = %lf",&ddum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->ordmedrej<0.0) par->ordmedrej=ddum;
  }
  READ_DATA_FILE;
  READ_DATA_FILE;
  if (par->clipsig<0.0 && sscanf(buffer," clipsig = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->clipsig<0.0) par->clipsig=ddum;
  READ_DATA_FILE;
  if (par->lrgerr<0.0 && sscanf(buffer," lrgerr = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->lrgerr<0.0) par->lrgerr=ddum;
  if (par->version>0.12) {
    READ_DATA_FILE;
    if (par->rankspec<0 && sscanf(buffer," rankspec = %d",&idum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->rankspec<0) par->rankspec=idum;
  }
  if (par->version>0.23) {
    READ_DATA_FILE;
    if (par->scalclip<0.0 && sscanf(buffer," scalmeth = %d",&idum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->scalmeth<0) par->scalmeth=idum;
  } else if (par->scalmeth<0) par->scalmeth=0;
  READ_DATA_FILE;
  if (par->scalclip<0.0 && sscanf(buffer," scalclip = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->scalclip<0.0) par->scalclip=ddum;
  READ_DATA_FILE;
  if (par->scalerr<0.0 && sscanf(buffer," scalerr = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->scalerr<0.0) par->scalerr=ddum;
  if (par->version>0.25) {
    READ_DATA_FILE;
    if (par->nscalclip<0 && sscanf(buffer," nscalclip = %d",&idum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->nscalclip<0) par->nscalclip=idum;
  } else if (par->nscalclip<0) par->nscalclip=NSCALCLIP;
  READ_DATA_FILE;
  READ_DATA_FILE;
  if (par->contftyp<0 && sscanf(buffer," contftyp = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->contftyp<0) par->contftyp=idum;
  READ_DATA_FILE;
  if (par->contord<0 && sscanf(buffer," contord = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->contord<0) par->contord=idum;
  READ_DATA_FILE;
  if (par->contpctl<0.0 && sscanf(buffer," contpctl = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->contpctl<0.0) par->contpctl=ddum;
  READ_DATA_FILE;
  if (par->contsigl<0.0 && sscanf(buffer," contsigl = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->contsigl<0.0) par->contsigl=ddum;
  READ_DATA_FILE;
  if (par->contsigu<0.0 && sscanf(buffer," contsigu = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->contsigu<0.0) par->contsigu=ddum;
  READ_DATA_FILE;
  if (par->contwgt<0 && sscanf(buffer," contwgt = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->contwgt<0) par->contwgt=idum;
  READ_DATA_FILE;
  READ_DATA_FILE;
  if (par->cordlya<0 && sscanf(buffer," cordlya = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->cordlya<0) par->cordlya=idum;
  READ_DATA_FILE;
  if (par->cordred<0 && sscanf(buffer," cordred = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->cordred<0) par->cordred=idum;
  READ_DATA_FILE;
  if (par->ftyplya<0 && sscanf(buffer," ftyplya = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->ftyplya<0) par->ftyplya=idum;
  READ_DATA_FILE;
  if (par->ftypred<0 && sscanf(buffer," ftypred = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->ftypred<0) par->ftypred=idum;
  if (par->version>0.15) {
    READ_DATA_FILE;
    if (par->nocont<0 && sscanf(buffer," nocont = %d",&idum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->nocont<0) par->nocont=idum;
  }
  READ_DATA_FILE;
  if (par->pctllya<0.0 && sscanf(buffer," pctllya = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->pctllya<0.0) par->pctllya=ddum;
  READ_DATA_FILE;
  if (par->pctlred<0.0 && sscanf(buffer," pctlred = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->pctlred<0.0) par->pctlred=ddum; 
  READ_DATA_FILE;
  if (par->rsiglyal<0.0 && sscanf(buffer," rsiglyal = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->rsiglyal<0.0) par->rsiglyal=ddum;
  READ_DATA_FILE;
  if (par->rsiglyau<0.0 && sscanf(buffer," rsiglyau = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->rsiglyau<0.0) par->rsiglyau=ddum;
  READ_DATA_FILE;
  if (par->rsigredl<0.0 && sscanf(buffer," rsigredl = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->rsigredl<0.0) par->rsigredl=ddum;
  READ_DATA_FILE;
  if (par->rsigredu<0.0 && sscanf(buffer," rsigredu = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->rsigredu<0.0) par->rsigredu=ddum;
  READ_DATA_FILE;
  if (par->vclya<0.0 && sscanf(buffer," vclya = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->vclya<0.0) par->vclya=ddum;
  READ_DATA_FILE;
  if (par->vcred<0.0 && sscanf(buffer," vcred = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->vcred<0.0) par->vcred=ddum;
  READ_DATA_FILE;
  if (par->vlya<0.0 && sscanf(buffer," vlya = %lf",&ddum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->vlya<0.0) par->vlya=ddum;
  READ_DATA_FILE;
  READ_DATA_FILE;
  if (par->dat<0 && sscanf(buffer," dat = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->dat<0) par->dat=idum;
  if (par->version>0.52) {
    READ_DATA_FILE;
    if (par->raw<0 && sscanf(buffer," raw = %d",&idum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->raw<0) par->raw=idum;
  }
  READ_DATA_FILE;
  if (par->save<0 && sscanf(buffer," save = %d",&idum)!=1) {
    fclose(data_file); ERR_FMT;
  } else if (par->save<0) par->save=idum;
  if (par->version>=0.56) {
    READ_DATA_FILE;
    if (sscanf(buffer," distort = %d",&idum)!=1) {
      fclose(data_file); ERR_FMT;
    } else if (par->distort<0) par->distort=idum;
    else par->distort=(!idum) ? 1 : 0;
  }

  /* Get number of spectra */
  READ_DATA_FILE;
  READ_DATA_FILE;
  if (sscanf(buffer,"Number_of_spectra = %d",nspec)!=1) {
    fclose(data_file); ERR_FMT;
  }
  
  /* Allocate memory for required number of spectra */
  if (!(*spec=(spectrum *)malloc((size_t)(*nspec*sizeof(spectrum)))))
    errormsg("UVES_rUPLfile(): Could not allocate memory for\n\
\tspectrum array of size %d",*nspec);

  /* If the user wants to use a Macmap file to map the path and file
     names to something else, read in Macmap file here */
  if (par->macmap==1) {
    sprintf(mmap.mmapfile,"%s",par->macmapfile);
    if (!UVES_rMacmap(&mmap)) {
      nferrormsg("UVES_rUPLfile(): Error reading Macmap file %s",mmap.mmapfile);
      return 0;
    }
  } else mmap.nmap=0;
  
  /* Get file paths, file names, velocity shifts, random distortion
     seeds and assign identity numbers */
  /* Treat special case of reading in an already-combined spectrum first */
  if (par->filetype==FTCOMB) {
    READ_DATA_FILE;
    if (sscanf(buffer,"%*d_PATH = %s",cspec->path)!=1) {
      fclose(data_file); ERR_FMT;
    }
    READ_DATA_FILE;
    if (sscanf(buffer," %*d_COMB = %s",cspec->abfile)!=1) {
      fclose(data_file); ERR_FMT;
    }
    sprintf(cspec->file,"%s%s",cspec->path,cspec->abfile);
    READ_DATA_FILE;
  } else {
    /* Otherwise, read in details for each separate spectrum */
    READ_DATA_FILE;
    if (sscanf(buffer,"%d_%*s",&j)!=1) { fclose(data_file); ERR_FMT; }
    for (i=0; i<*nspec; i++) {
      (*spec)[i].ftype=-2; (*spec)[i].vshift=0.0; (*spec)[i].distort_seed=0;
      while (sscanf(buffer,"%d_%*s",&k)==1 && k==j) {
	if (sscanf(buffer,"%*d_FTYP = %d",&((*spec)[i].ftype))==1) READ_DATA_FILE;
	if (sscanf(buffer,"%*d_PATH = %s",(*spec)[i].path)==1) READ_DATA_FILE;
	if (sscanf(buffer," %*d_FLUX = %s",(*spec)[i].abfile)==1) READ_DATA_FILE;
	if (sscanf(buffer," %*d_ERRO = %s",(*spec)[i].aberfile)==1) READ_DATA_FILE;
	if (sscanf(buffer," %*d_THAR = %s",(*spec)[i].abthfile)==1) READ_DATA_FILE;
	if (sscanf(buffer," %*d_WPOL = %s",(*spec)[i].abwlfile)==1) READ_DATA_FILE;
	if (sscanf(buffer," %*d_VSHT = %lf",&((*spec)[i].vshift))==1) READ_DATA_FILE;
	if (sscanf(buffer," %*d_DTRT = %ld",&((*spec)[i].distort_seed))==1)
	  READ_DATA_FILE;
      }
      j=k;
    }
  }
  /* Check file paths, file names etc. for each spectrum to ensure
     they're mutually consistent */
  for (i=0; i<*nspec; i++) {
    /* Ensure filetype is set */
    if (par->filetype==FTMIX && (*spec)[i].ftype==-2) {
      fclose(data_file);
      errormsg("UVES_rUPLfile(): No filetype flag read from UPL file\n\
\t%s for spectrum %d\n\t%s.\n\tMixed filetype was specified so each\n\
\tspectrum must have individual filetype flag in UPL file.",infile,i+1,
	       (*spec)[i].file);
    } else if (par->filetype==-2) (*spec)[i].ftype=FILETYPE;
    else if (par->filetype!=FTMIX && (*spec)[i].ftype==-2)
      (*spec)[i].ftype=par->filetype;
    /* Ensure path and necessary files have been specified */
    if (strlen((*spec)[i].path)<3)
      errormsg("UVES_rUPLfile(): No path read from UPL file\n\
\t%s for spectrum %d\n\t%s",infile,i+1,(*spec)[i].file);
    if (strlen((*spec)[i].abfile)<3)
      errormsg("UVES_rUPLfile(): No flux file read from UPL file\n\
\t%s for spectrum %d\n\t%s",infile,i+1,(*spec)[i].file);
    if (par->thar<=1) {
      if (strlen((*spec)[i].aberfile)<3)
	errormsg("UVES_rUPLfile(): No flux file read from UPL file\n\
\t%s for spectrum %d\n\t%s",infile,i+1,(*spec)[i].file);
    }
    if (strlen((*spec)[i].abwlfile)<3)
      errormsg("UVES_rUPLfile(): No wavelength info file read from UPL file\n\
\t%s for spectrum %d\n\t%s",infile,i+1,(*spec)[i].file);
    if (par->thar==1) {
      /* Set default file names for the ThAr spectra if the user has
	 not entered them in the UPL file */
      if ((cptr=strstr((*spec)[i].abfile,"_"))==NULL)
	errormsg("UVES_rUPLfile(): Cannot identify an '_' in flux file\n\
\tname\n\t%s\n\tto separate file type and object name. I therefore cannot\n \
\tdetermine the default name for the associated ThAr file. Please rename your\n	\
\tfiles according to the correct conventions and try again.",(*spec)[i].abfile);
      cptr++;
      switch ((*spec)[i].ftype) {
      case FTUVES:
	sprintf((*spec)[i].abthfile,"thar_%s",cptr); break;
      case FTIRAF:
	sprintf((*spec)[i].abthfile,"wav_%s",cptr); break;
      case FTMAKE:
	sprintf((*spec)[i].abthfile,"Wav-%s",cptr); break;
      case FTIRLS:
	sprintf((*spec)[i].abthfile,"wav_%s",cptr); break;
      case FTHIRX:
	sprintf((*spec)[i].abthfile,"arc_%s",cptr); break;
      }
    }
    /** If the user wants to use a Macmap file to map the path and file
       names to something else, do the mapping here **/
    /* Match flux file name to Macmap parameters */
    for (j=0; j<mmap.nmap; j++) {
      sprintf(stro,"%s_%02d",mmap.cwl_ori[j],mmap.idx_ori[j]);
      sprintf(newname,"%s",(*spec)[i].abfile);
      if ((cp1=strstr(newname,mmap.obj_ori[j])) && !strncmp(cp1-1,"_",1) &&
	  !strncmp((cp2=cp1+strlen(mmap.obj_ori[j])),"_",1) && 
	  (cp3=strstr(newname,stro))) {
	/* Map the flux file name */
	cp4=cp3+strlen(stro); sprintf(strn,"%s_%02d",mmap.cwl_new[j],mmap.idx_new[j]);
	strncpy((*spec)[i].abfile,newname,cp1-newname);
	(*spec)[i].abfile[cp1-newname]='\0'; strcat((*spec)[i].abfile,mmap.obj_new[j]);
	strncat((*spec)[i].abfile,cp2,cp3-cp2); strcat((*spec)[i].abfile,strn);
	strcat((*spec)[i].abfile,cp4); strcat((*spec)[i].abfile,"\0");
	/* Map the error file name */
	sprintf(newname,"%s",(*spec)[i].aberfile); cp1=strstr(newname,mmap.obj_ori[j]);
	cp2=cp1+strlen(mmap.obj_ori[j]); cp3=strstr(newname,stro); cp4=cp3+strlen(stro);
	strncpy((*spec)[i].aberfile,newname,cp1-newname);
	(*spec)[i].aberfile[cp1-newname]='\0'; strcat((*spec)[i].aberfile,mmap.obj_new[j]);
	strncat((*spec)[i].aberfile,cp2,cp3-cp2); strcat((*spec)[i].aberfile,strn);
	strcat((*spec)[i].aberfile,cp4); strcat((*spec)[i].aberfile,"\0");
	/* Map the wavelength polynomial file name */
	sprintf(newname,"%s",(*spec)[i].abwlfile); cp1=strstr(newname,mmap.obj_ori[j]);
	cp2=cp1+strlen(mmap.obj_ori[j]); cp3=strstr(newname,stro); cp4=cp3+strlen(stro);
	strncpy((*spec)[i].abwlfile,newname,cp1-newname);
	(*spec)[i].abwlfile[cp1-newname]='\0'; strcat((*spec)[i].abwlfile,mmap.obj_new[j]);
	strncat((*spec)[i].abwlfile,cp2,cp3-cp2); strcat((*spec)[i].abwlfile,strn);
	strcat((*spec)[i].abwlfile,cp4); strcat((*spec)[i].abwlfile,"\0");
	/* Map the ThAr file name if necessary */
	if (par->thar==1) {
	  sprintf(newname,"%s",(*spec)[i].abthfile); cp1=strstr(newname,mmap.obj_ori[j]);
	  cp2=cp1+strlen(mmap.obj_ori[j]); cp3=strstr(newname,stro); cp4=cp3+strlen(stro);
	  strncpy((*spec)[i].abthfile,newname,cp1-newname);
	  (*spec)[i].abthfile[cp1-newname]='\0'; strcat((*spec)[i].abthfile,mmap.obj_new[j]);
	  strncat((*spec)[i].abthfile,cp2,cp3-cp2); strcat((*spec)[i].abthfile,strn);
	  strcat((*spec)[i].abthfile,cp4); strcat((*spec)[i].abthfile,"\0");
	}
	/* Map the path name */
	sprintf(newname,"%s",(*spec)[i].path); cp1=NULL; cp2=newname;
	while ((cp2=strstr(cp2,mmap.obj_ori[j]))!=NULL) { cp1=cp2; cp2++; }
	cp2=(cp1) ? cp1+strlen(mmap.obj_ori[j]) : NULL;
	if (cp1 && !strncmp(cp1-1,"/",1) && !strncmp(cp2,"/",1)) {
	  strncpy((*spec)[i].path,newname,cp1-newname); (*spec)[i].path[cp1-newname]='\0';
	  strcat((*spec)[i].path,mmap.obj_new[j]); strcat((*spec)[i].path,cp2);
	  strcat((*spec)[i].path,"\0");
	}
	break;
      }
    }
    /* Construct full file names by concatenating path and abbreviated
       file names */
    sprintf((*spec)[i].file,"%s%s",(*spec)[i].path,(*spec)[i].abfile);
    sprintf((*spec)[i].erfile,"%s%s",(*spec)[i].path,(*spec)[i].aberfile);
    sprintf((*spec)[i].wlfile,"%s%s",(*spec)[i].path,(*spec)[i].abwlfile);
    if (par->thar==1)
      sprintf((*spec)[i].thfile,"%s%s",(*spec)[i].path,(*spec)[i].abthfile);
    /* Assign spectrum identity numbers */
    (*spec)[i].id=i;
  }
  
  /* Get number of actions */
  READ_DATA_FILE;
  if (sscanf(buffer,"Number_of_actions = %d",nact)!=1) {
    fclose(data_file); ERR_FMT;
  }
  
  /* Allocate memory for required number of actions */
  if (!(*act=(action *)malloc((size_t)(*nact*sizeof(action)))))
    errormsg("UVES_rUPLfile(): Could not allocate memory for\n\
\taction array of size %d",*nact);
  
  /* Get action details */
  for (i=0; i<*nact; i++) {
    READ_DATA_FILE;
    if (sscanf(buffer,"%*d_ACTN = %d  %*d_RCMB = %d",&((*act)[i].act),
	       &((*act)[i].rcmb))!=2) {
      fclose(data_file); ERR_FMT;
    }
    if ((*act)[i].act<=ICACT) {
      READ_DATA_FILE;
      if (sscanf(buffer," %*d_CORD = %lf %lf %lf %lf",&((*act)[i].d[0]),
		 &((*act)[i].d[1]),&((*act)[i].d[2]),&((*act)[i].d[3]))!=4) {
	fclose(data_file); ERR_FMT;
      }
    }
    
    /* Some action-specific entries */
    switch((*act)[i].act) {
    case COACT:
      READ_DATA_FILE;
      if (sscanf(buffer," %*d_SPEC = %d  %*d_ORDR = %d",&((*act)[i].i[0]),
		 &((*act)[i].i[1]))!=2) {
	fclose(data_file); ERR_FMT;
      }
      (*act)[i].i[0]--; (*act)[i].i[1]--;
      break;
    case UOACT:
      READ_DATA_FILE;
      if (sscanf(buffer," %*d_SPEC = %d  %*d_ORDR = %d",&((*act)[i].i[0]),
		 &((*act)[i].i[1]))!=2) {
	fclose(data_file); ERR_FMT;
      }
      (*act)[i].i[0]--; (*act)[i].i[1]--;
      break;
    case FOACT:
      READ_DATA_FILE;
      if (sscanf(buffer," %*d_SPEC = %d  %*d_ORDR = %d  %*d_FTYP = %d  \
%*d_FORD = %d  %*d_REJL = %lf  %*d_REJU = %lf  %*d_NXYP = %d",&((*act)[i].i[0]),
		 &((*act)[i].i[1]),&((*act)[i].i[2]),&((*act)[i].i[3]),
		 &((*act)[i].d[4]),&((*act)[i].d[5]),&((*act)[i].nxyp))!=7) {
	fclose(data_file); ERR_FMT;
      }
      (*act)[i].i[0]--; (*act)[i].i[1]--;
      if ((*act)[i].nxyp) {
	/* Allocate memory for any xy-pairs */
	if (!((*act)[i].xyp=
	      (twoxyp *)malloc((size_t)((*act)[i].nxyp*sizeof(twoxyp))))) {
	  fclose(data_file);
	  errormsg("UVES_rUPLfile(): Cannot allocate memory for\n\
\txyp array of size %d",(*act)[i].nxyp);
	}
	/* Get any xy-pairs */
	for (j=0; j<(*act)[i].nxyp; j++) {
	  READ_DATA_FILE;
	  if (sscanf(buffer,"  %*d_XY_%*d = %d %lf %lf %lf %lf",
		     &((*act)[i].xyp[j].i),&((*act)[i].xyp[j].x1),
		     &((*act)[i].xyp[j].x2),&((*act)[i].xyp[j].y1),
		     &((*act)[i].xyp[j].y2))!=5) {
	    fclose(data_file); ERR_FMT;
	  }
	}
      }
      break;
    case FCACT:
      READ_DATA_FILE;
      if (sscanf(buffer," %*d_FTYP = %d  %*d_FORD = %d  %*d_REJL = %lf  \
%*d_REJU = %lf  %*d_NXYP = %d",&((*act)[i].i[0]),&((*act)[i].i[1]),
		 &((*act)[i].d[4]),&((*act)[i].d[5]),&((*act)[i].nxyp))!=5) {
	fclose(data_file); ERR_FMT;
      }
      if ((*act)[i].nxyp) {
	/* Allocate memory for any xy-pairs */
	if (!((*act)[i].xyp=
	      (twoxyp *)malloc((size_t)((*act)[i].nxyp*sizeof(twoxyp))))) {
	  fclose(data_file);
	  errormsg("UVES_rUPLfile(): Cannot allocate memory for\n\
\txyp array of size %d",(*act)[i].nxyp);
	}
	/* Get any xy-pairs */
	for (j=0; j<(*act)[i].nxyp; j++) {
	  READ_DATA_FILE;
	  if (sscanf(buffer,"  %*d_XY_%*d = %d %lf %lf %lf %lf",
		     &((*act)[i].xyp[j].i),&((*act)[i].xyp[j].x1),
		     &((*act)[i].xyp[j].x2),&((*act)[i].xyp[j].y1),
		     &((*act)[i].xyp[j].y2))!=5) {
	    fclose(data_file); ERR_FMT;
	  }
	}
      }
      break;
    case IOACT:
      READ_DATA_FILE;
      if (sscanf(buffer," %*d_SPEC = %d  %*d_ORDR = %d  %*d_NXYP = %d",
		 &((*act)[i].i[0]),&((*act)[i].i[1]),&((*act)[i].nxyp))!=3) {
	fclose(data_file); ERR_FMT;
      }
      (*act)[i].i[0]--; (*act)[i].i[1]--;
      if (!((*act)[i].nxyp)) {
	fclose(data_file);
	errormsg("UVES_rUPLfile(): Error reading line %d in file %s.\n\
\tNumber of XY pairs must be >=2",idx,infile);
      }
      /* Allocate memory for any xy-pairs */
      if (!((*act)[i].xyp=
	    (twoxyp *)malloc((size_t)((*act)[i].nxyp*sizeof(twoxyp))))) {
	fclose(data_file);
	errormsg("UVES_rUPLfile(): Cannot allocate memory for\n\
\txyp array of size %d",(*act)[i].nxyp);
      }
      /* Get any xy-pairs */
      for (j=0; j<(*act)[i].nxyp; j++) {
	READ_DATA_FILE;
	if (sscanf(buffer,"  %*d_XY_%*d = %lf %lf",&((*act)[i].xyp[j].x1),
		   &((*act)[i].xyp[j].y1))!=2) {
	  fclose(data_file); ERR_FMT;
	}
      }
      break;
    case ICACT:
      READ_DATA_FILE;
      if (sscanf(buffer," %*d_NXYP = %d",&((*act)[i].nxyp))!=1) {
	fclose(data_file); ERR_FMT;
      }
      if (!((*act)[i].nxyp)) {
	fclose(data_file);
	errormsg("UVES_rUPLfile(): Error reading line %d in file %s.\n\
\tNumber of XY pairs must be >=2",idx,infile);
      }
      /* Allocate memory for any xy-pairs */
      if (!((*act)[i].xyp=
	    (twoxyp *)malloc((size_t)((*act)[i].nxyp*sizeof(twoxyp))))) {
	fclose(data_file);
	errormsg("UVES_rUPLfile(): Cannot allocate memory for\n\
\txyp array of size %d",(*act)[i].nxyp);
      }
      /* Get any xy-pairs */
      for (j=0; j<(*act)[i].nxyp; j++) {
	READ_DATA_FILE;
	if (sscanf(buffer,"  %*d_XY_%*d = %lf %lf",&((*act)[i].xyp[j].x1),
		   &((*act)[i].xyp[j].y1))!=2) {
	  fclose(data_file); ERR_FMT;
	}
      }
      break;
    case SOACT:
      READ_DATA_FILE;
      if (sscanf(buffer," %*d_SPEC = %d  %*d_ORDR = %d  %*d_SCAL = %lf",
		 &((*act)[i].i[0]),&((*act)[i].i[1]),&((*act)[i].d[0]))!=3) {
	fclose(data_file); ERR_FMT;
      }
      (*act)[i].i[0]--; (*act)[i].i[1]--;
      break;
    case ARACT:
      READ_DATA_FILE;
      if (sscanf(buffer," %*d_SCLP = %lf  %*d_SERR = %lf",
		 &((*act)[i].d[0]),&((*act)[i].d[1]))!=2) {
	fclose(data_file); ERR_FMT;
      }
      break;
    }
    /* Make every action a valid one */
    (*act)[i].val=1;
  }

  /* If we're not at the end of the file here then there might be
     actions that haven't been read in */
  cptr=fgets(buffer,LNGSTRLEN,data_file);
  if (!feof(data_file))
    warnmsg("UVES_rUPLfile(): Read in %d actions but more actions\n\
\tmay exist in UPL file, %s. I am continuing but, be warned,\n		\
\tyou may lose actions. Press ctrl-c to quit now if desired.",*nact,infile);
  fclose(data_file);

  return 1;
  
}
