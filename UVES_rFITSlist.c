/****************************************************************************
* Read in list of FITS files with absolute path names
****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "UVES_popler.h"
#include "file.h"
#include "error.h"

int UVES_rFITSlist(char *infile, spectrum **spec, int *nspec, params *par) {

  int    i=0;
  long   ranseed=0;
  char   buffer[LNGSTRLEN]="\0";
  char   *cptr;
  FILE   *data_file=NULL;

  /* Open input file */
  if ((data_file=faskropen("Valid input UPL file, FITS file, FITS file list?",
			   infile,5))==NULL)
    errormsg("UVES_rFITSlist(): Can not open file %s",infile);
  /* Read in first line  */
  i++; if ((cptr=fgets(buffer,LNGSTRLEN,data_file))==NULL) {
    fclose(data_file);
    errormsg("UVES_rFITSlist(): Problem reading file %s on line %d",infile,i);
  }
  /* File should be a list of FITS file names with absolute paths */
  /* Check for absolute path names or environment variables ... a weak check anyway */
  if (strncmp(buffer,"/",1) && strncmp(buffer,"$",1)) {
    fclose(data_file);
      errormsg("UVES_rFITSlist(): FITS file path invalid on line %d\n\
\tin file %s.\n\
\tYou must use absolute path names directly or in environment\n\
\tvariable names, i.e. FITS file names must begin with '/' or '$'.",
	       1,infile);
  }
  /* See if it is a list of valid FITS files */
  while ((cptr=fgets(buffer,LNGSTRLEN,data_file))!=NULL) {
    /* Check for absolute path names or environment variables ... a weak check anyway */
    if (strncmp(buffer,"/",1) && strncmp(buffer,"$",1)) {
      fclose(data_file);
      errormsg("UVES_rFITSlist(): FITS file path invalid on line %d\n\
\tin file %s.\n\
\tYou must use absolute path names - FITS file names must begin with '/'.",
	       i+1,infile);
    }
    i++;
  }
  if (!feof(data_file)) {
    fclose(data_file);
    errormsg("UVES_rFITSlist(): Error reading line %d in file %s",i,infile);
  }
  else {
    rewind(data_file); *nspec=i;
    /* Allocate memory for nspec spectra */
    if (!(*spec=(spectrum *)malloc((size_t)(*nspec*sizeof(spectrum)))))
      errormsg("UVES_rFITSlist(): Could not allocate memory for\n\
\tspectrum array of size %d",*nspec);
  }

  /* Read in list of names of FITS files and assign identity numbers */
  for (i=0; i<*nspec; i++) {
    cptr=fgets(buffer,LNGSTRLEN,data_file);
    if (par->filetype==FTMIX) {
      if (sscanf(buffer,"%s %d",(*spec)[i].file,&((*spec)[i].ftype))!=2) {
	fclose(data_file);
	errormsg("UVES_rFITSlist(): Incorrect format in line %d\n\t of file\n\
\t%s.\n\tWhen using mixed file-types, you must enter a digit (%d-%d) after each\n\
\tfile name of the input file",i+1,infile,FTUVES,FTHARP);
      }
    } else if (sscanf(buffer,"%s",(*spec)[i].file)!=1) {
      fclose(data_file);
      errormsg("UVES_rFITSlist(): Incorrect format in line %d\n\t of file %s",
	       i+1,infile);
    } else if (par->filetype==-2) (*spec)[i].ftype=FILETYPE;
    else (*spec)[i].ftype=par->filetype;
    /* Get file name and path etc. */
    cptr=strrchr((*spec)[i].file,'/'); strcpy((*spec)[i].abfile,cptr+1);
    cptr=strrchr(strcpy((*spec)[i].path,(*spec)[i].file),'/'); *(cptr+1)='\0';
    /* Assign identity number */
    (*spec)[i].id=i;
  }
  fclose(data_file);

  /* Generate list of error and wpol file names */
  for (i=0; i<*nspec; i++) {
    if (par->thar<=1) {
      if ((*spec)[i].ftype<=FTUVES) {
	if ((cptr=strstr((*spec)[i].abfile,"fxb_"))==NULL ||
	    strncmp((*spec)[i].abfile,"fxb_",4))
	  errormsg("UVES_rFITSlist(): Cannot identify error and\n\
\twavelength polynomial array files for flux file\n\t%s\n\
\tUVES flux files must have 'fxb_' prefix",(*spec)[i].file);
	cptr+=4;
	sprintf((*spec)[i].erfile,"%serr_%s",(*spec)[i].path,cptr);
	sprintf((*spec)[i].aberfile,"err_%s",cptr);
	if (par->thar==1) {
	  sprintf((*spec)[i].thfile,"%sthar_%s",(*spec)[i].path,cptr);
	  sprintf((*spec)[i].abthfile,"thar_%s",cptr);
	}
	sprintf((*spec)[i].wlfile,"%swpol_%s",(*spec)[i].path,cptr);
	sprintf((*spec)[i].abwlfile,"wpol_%s",cptr);
      } else if ((*spec)[i].ftype==FTIRAF) {
	if ((cptr=strstr((*spec)[i].abfile,"flx_"))==NULL ||
	    strncmp((*spec)[i].abfile,"flx_",4))
	  errormsg("UVES_rFITSlist(): Cannot identify error and\n\
\twavelength polynomial array files for flux file\n\t%s\n\
\tIRAF flux files must have 'flx_' prefix",(*spec)[i].file);
	cptr+=4;
	sprintf((*spec)[i].erfile,"%ssig_%s",(*spec)[i].path,cptr);
	sprintf((*spec)[i].aberfile,"sig_%s",cptr);
	if (par->thar==1) {
	  sprintf((*spec)[i].thfile,"%swav_%s",(*spec)[i].path,cptr);
	  sprintf((*spec)[i].abthfile,"wav_%s",cptr);
	}
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTMAKE) {
	if ((cptr=strstr((*spec)[i].abfile,"Flux-"))==NULL ||
	    strncmp((*spec)[i].abfile,"Flux-",5))
	  errormsg("UVES_rFITSlist(): Cannot identify error and\n\
\twavelength polynomial array files for flux file\n\t%s\n\
\tMAKEE flux files must have 'Flux-' prefix",(*spec)[i].file);
	cptr+=5;
	sprintf((*spec)[i].erfile,"%sErr-%s",(*spec)[i].path,cptr);
	sprintf((*spec)[i].aberfile,"Err-%s",cptr);
	if (par->thar==1) {
	  sprintf((*spec)[i].thfile,"%sWav-%s",(*spec)[i].path,cptr);
	  sprintf((*spec)[i].abthfile,"Wav-%s",cptr);
	}
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTIRLS) {
	if ((cptr=strstr((*spec)[i].abfile,"flx_"))==NULL ||
	    strncmp((*spec)[i].abfile,"flx_",4))
	  errormsg("UVES_rFITSlist(): Cannot identify error and\n\
\twavelength polynomial array files for flux file\n\t%s\n\
\tIRAFLSS flux files must have 'flx_' prefix",(*spec)[i].file);
	cptr+=4;
	sprintf((*spec)[i].erfile,"%ssig_%s",(*spec)[i].path,cptr);
	sprintf((*spec)[i].aberfile,"sig_%s",cptr);
	if (par->thar==1) {
	  sprintf((*spec)[i].thfile,"%swav_%s",(*spec)[i].path,cptr);
	  sprintf((*spec)[i].abthfile,"wav_%s",cptr);
	}
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTHIRX) {
	sprintf((*spec)[i].erfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].aberfile,"%s",(*spec)[i].abfile);
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
	if (par->thar==1) {
	  if ((cptr=strstr((*spec)[i].abfile,"spec_"))==NULL ||
	      strncmp((*spec)[i].abfile,"spec_",5))
	    errormsg("UVES_rFITSlist(): Cannot identify extracted ThAr\n\
\tarray files for flux file\n\t%s\n\tHIREDUX flux files must have 'spec_' prefix",
		     (*spec)[i].file);
	  cptr+=5;
	  sprintf((*spec)[i].thfile,"%sarc_%s",(*spec)[i].path,cptr);
	  sprintf((*spec)[i].abthfile,"arc_%s",cptr);
	}
      } else if ((*spec)[i].ftype==FTESOM) {
	if ((cptr=strstr((*spec)[i].abfile,"UV_SRED_"))==NULL ||
	    strncmp((*spec)[i].abfile,"UV_SRED_",8))
	  errormsg("UVES_rFITSlist(): Cannot identify error and\n\
\twavelength polynomial array files for flux file\n\t%s\n\
\tUVES ESO-merged flux files must have 'UV_SRED_' prefix",(*spec)[i].file);
	cptr+=8;
	sprintf((*spec)[i].erfile,"%sUV_SERR_%s",(*spec)[i].path,cptr);
	sprintf((*spec)[i].aberfile,"UV_SERR_%s",cptr);
	if (par->thar==1) {
	  sprintf((*spec)[i].thfile,"%sthar_%s",(*spec)[i].path,cptr);
	  sprintf((*spec)[i].abthfile,"thar_%s",cptr);
	}
	sprintf((*spec)[i].wlfile,"%sUV_SRED_%s",(*spec)[i].path,cptr);
	sprintf((*spec)[i].abwlfile,"UV_SRED_%s",cptr);
      } else if ((*spec)[i].ftype==FTKODI) {
	/* At the moment, KODIAQ files can have any file name because
	   no associations must be assumed between flux files and
	   other files (e.g. error, wavelength); the flux, wavelength
	   and error information is all in one file */
	sprintf((*spec)[i].erfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].aberfile,"%s",(*spec)[i].abfile);
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTPYPE) {
	/* At the moment, PypeIt files can have any file name because
	   no associations must be assumed between flux files and
	   other files (e.g. error, wavelength); the flux, wavelength
	   and error information is all in one file */
	sprintf((*spec)[i].erfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].aberfile,"%s",(*spec)[i].abfile);
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTMAGE) {
	sprintf((*spec)[i].aberfile,"%s",(*spec)[i].abfile);
	if ((cptr=strstr((*spec)[i].aberfile,"spec_cal"))==NULL)
	  errormsg("UVES_rFITSlist(): Cannot identify error and\n\
\twavelength polynomial array files for flux file\n\t%s\n\
\tMagE flux files must have 'spec_cal' in their names",(*spec)[i].file);
	sprintf(buffer,"%s",cptr+8); sprintf(cptr,"%s","err_cal");
	sprintf(cptr+7,"%s",buffer);
	sprintf((*spec)[i].erfile,"%s%s",(*spec)[i].path,(*spec)[i].aberfile);
	if (par->thar==1) {
	  sprintf((*spec)[i].abthfile,"%s",(*spec)[i].abfile);
	  cptr=strstr((*spec)[i].abthfile,"spec_cal");
	  sprintf(cptr,"%s","thar"); sprintf(cptr+4,"%s",buffer);
	  sprintf((*spec)[i].thfile,"%s%s",(*spec)[i].path,(*spec)[i].abthfile);
	}
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
	cptr=strstr((*spec)[i].abwlfile,"spec_cal");
	sprintf(cptr,"%s","wspec"); sprintf(cptr+5,"%s",buffer);
	sprintf((*spec)[i].wlfile,"%s%s",(*spec)[i].path,(*spec)[i].abwlfile);
      } else if ((*spec)[i].ftype==FTIESI) {
	sprintf((*spec)[i].erfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].aberfile,"%s",(*spec)[i].abfile);
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTHARP) {
	if ((cptr=strstr((*spec)[i].abfile,"HARP"))==NULL ||
	    strncmp((*spec)[i].abfile,"HARP",4) ||
	    (cptr=strstr((*spec)[i].abfile,"e2ds"))==NULL)
	  errormsg("UVES_rFITSlist(): Cannot identify HARPS filename\n\
\tcharacteristics for flux file\n\t%s\n\
\tHARPS (S or N) file names must have 'HARP[S or N]' prefix\n\
\tand contain 'e2ds'.",(*spec)[i].file);
	sprintf((*spec)[i].erfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].aberfile,"%s",(*spec)[i].abfile);
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTESPR) {
	sprintf((*spec)[i].erfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].aberfile,"%s",(*spec)[i].abfile);
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      }
    } else {
      if ((*spec)[i].ftype<=FTUVES) {
	if ((cptr=strstr((*spec)[i].abfile,"thar_"))==NULL ||
	    strncmp((*spec)[i].abfile,"thar_",5))
	  errormsg("UVES_rFITSlist(): Cannot identify extracted\n\
\tThAr files for flux file\n\t%s\n\
\tUVES ThAr file must have 'thar_' prefix",(*spec)[i].file);
	cptr+=5;
	sprintf((*spec)[i].wlfile,"%swpol_%s",(*spec)[i].path,cptr);
	sprintf((*spec)[i].abwlfile,"wpol_%s",cptr);
      } else if ((*spec)[i].ftype==FTIRAF) {
	if ((cptr=strstr((*spec)[i].abfile,"wav_"))==NULL ||
	    strncmp((*spec)[i].abfile,"wav_",4))
	  errormsg("UVES_rFITSlist(): Cannot identify extracted\n\
\tThAr files for flux file\n\t%s\n\
\tIRAF ThAr file must have 'wav_' prefix",(*spec)[i].file);
	cptr+=4;
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTMAKE) {
	if ((cptr=strstr((*spec)[i].abfile,"Wav-"))==NULL ||
	    strncmp((*spec)[i].abfile,"Wav-",4))
	  errormsg("UVES_rFITSlist(): Cannot identify extracted\n\
\tThAr flux files for flux file\n\t%s\n\
\tMAKEE ThAr file must have 'Wav-' prefix",(*spec)[i].file);
	cptr+=4;
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTIRLS) {
	if ((cptr=strstr((*spec)[i].abfile,"wav_"))==NULL ||
	    strncmp((*spec)[i].abfile,"wav_",4))
	  errormsg("UVES_rFITSlist(): Cannot identify extracted\n\
\tThAr flux files for flux file\n\t%s\n\
\tIRAFLSS ThAr file must have 'wav_' prefix",(*spec)[i].file);
	cptr+=4;
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTHIRX) {
	if ((cptr=strstr((*spec)[i].abfile,"Arc_"))==NULL ||
	    strncmp((*spec)[i].abfile,"Arc_",4))
	  errormsg("UVES_rFITSlist(): Cannot identify extracted\n\
\tThAr flux files for flux file\n\t%s\n\
\tHIREDUX ThAr file must have 'Arc_' prefix",(*spec)[i].file);
	cptr+=4;
	sprintf((*spec)[i].wlfile,"%s",(*spec)[i].file);
	sprintf((*spec)[i].abwlfile,"%s",(*spec)[i].abfile);
      } else if ((*spec)[i].ftype==FTESOM) {
	if ((cptr=strstr((*spec)[i].abfile,"thar_"))==NULL ||
	    strncmp((*spec)[i].abfile,"thar_",5))
	  errormsg("UVES_rFITSlist(): Cannot identify extracted\n\
\tThAr files for flux file\n\t%s\n\
\tUVES ThAr file must have 'thar_' prefix",(*spec)[i].file);
	cptr+=5;
	sprintf((*spec)[i].wlfile,"%swpol_%s",(*spec)[i].path,cptr);
	sprintf((*spec)[i].abwlfile,"wpol_%s",cptr);
      } else if ((*spec)[i].ftype==FTKODI) {
	errormsg("UVES_rFITSlist(): Do not know how to read in\n\
\tThAr KODIAQ files for flux file\n\t%s",(*spec)[i].file);
      } else if ((*spec)[i].ftype==FTPYPE) {
	errormsg("UVES_rFITSlist(): Do not know how to read in\n\
\tThAr PypeIt files for flux file\n\t%s",(*spec)[i].file);
      } else if ((*spec)[i].ftype==FTMAGE) {
	errormsg("UVES_rFITSlist(): Do not know how to read in\n\
\tThAr MagE files for flux file\n\t%s",(*spec)[i].file);
      } else if ((*spec)[i].ftype==FTIESI) {
	errormsg("UVES_rFITSlist(): Do not know how to read in\n\
\tThAr IRAF-ESI files for flux file\n\t%s",(*spec)[i].file);
      } else if ((*spec)[i].ftype==FTHARP) {
	errormsg("UVES_rFITSlist(): Do not know how to read in\n\
\tThAr HARPS spectrum for flux file\n\t%s",(*spec)[i].file);
      }
    }
  }

  /* Initialize the user-supplied velocity shift and scaling
     information to be applied to each spectrum */
  for (i=0; i<*nspec; i++) {
    (*spec)[i].vshift=(*spec)[i].vslope=(*spec)[i].refwav=0.0;
    strcpy((*spec)[i].inscl,"\0");
  }

  /* Initialize the random number seeds for the wavelength distortions */
  ranseed=-time(NULL);
  for (i=0; i<*nspec; i++) (*spec)[i].distort_seed=ranseed++;

  return 1;

}
