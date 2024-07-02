/****************************************************************************
* Read in the input file specified by the user and decide whether it's a
* UPL file or a list of FITS files
****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "UVES_popler.h"
#include "file.h"
#include "error.h"

int UVES_rinputfile(char *infile, spectrum **spec, int *nspec, action **act,
		    int *nact, cspectrum *cspec, params *par) {

  int    idx=0;
  int    i=0;
  char   *cptr;
  char   buffer[VLNGSTRLEN]="\0";
  FILE   *data_file=NULL;

  /** Open input file, see if it's a UPL file or a list of FITS files
      or a single FITS file **/
  if ((data_file=faskropen("Valid input UPL file, FITS file, FITS file list?",
			   infile,5))==NULL)
    errormsg("UVES_rinputfile(): Can not open file %s",infile);
  
  /* Read in first line from input file to decide its nature */
  idx++; if (fgets(buffer,VLNGSTRLEN,data_file)==NULL) {
    fclose(data_file);
    errormsg("UVES_rinputfile(): Problem reading file %s on line %d",infile,idx);
  }
  if (!strncmp(buffer,"SIMPLE  =",8) || strisnum(buffer,VLNGSTRLEN,0)) {
    /* Single file specified that looks suspiciously like a FITS file
       or an ASCII data file */
    fclose(data_file);
    /* No contributing spectra are to be used, only a combined spectrum */
    *nspec=0; par->filetype=FTCOMB;
    /* Initialize UPL filename */
    strcpy(cspec->UPLfile,"\0");
    /* Set the combined spectrum file name and path */
    sprintf(cspec->file,"%s",infile);
    if ((cptr=strrchr(infile,'/'))==NULL) {
      if (getcwd(cspec->path,VLNGSTRLEN)==NULL)
	errormsg("UVES_rinputfile(): Could not get path of\n\
\tcurrent working directory");
      sprintf(cspec->abfile,"%s",infile);
      sprintf(cspec->file,"%s/%s",cspec->path,infile); strcat(cspec->path,"/");
    } else {
      if (strncmp(cspec->file,"/",1))
	errormsg("UVES_rinputfile(): Input file must be in\n\
\tcurrent working directory or you must specify the absolute\n\
\tpathname of the file");
      sprintf(cspec->abfile,"%s",cptr+1); sprintf(cspec->path,"%s",infile);
      cptr=strrchr(cspec->path,'/'); strcpy(cptr+1,"\0");
    }
  } else if (!strncmp(buffer,"UVES_popler : Version ",22)) {
    /* File looks suspiciously like a UPL file */
    fclose(data_file);
    /* Read the UPL file */
    if (!UVES_rUPLfile(infile,spec,nspec,act,nact,cspec,par))
      errormsg("UVES_rinputfile(): Error reading suspected UPL file\n\t%s",infile);
    /* Record UPL filename */
    sprintf(cspec->UPLfile,"%s",infile);
  } else {
    /* File should be a list of FITS files names with absolute paths */
    fclose(data_file);
    /* Read the FITS file list */
    if (!UVES_rFITSlist(infile,spec,nspec,par))
      errormsg("UVES_rinputfile(): Error reading FITS file list\n\t%s",infile);
    /* Initialize UPL filename */
    strcpy(cspec->UPLfile,"\0");
  }

  /* Check for file existence */
  if (par->filetype==FTCOMB) {
    if (access(cspec->file,R_OK))
      errormsg("UVES_rinputfile(): Cannot find previously combined\n\
\tspectrum file \n\t%s",cspec->file);
  }
  for (i=0; i<*nspec; i++) {
    if (access(UVES_replace_envinstr((*spec)[i].file),R_OK))
      errormsg("UVES_rinputfile(): Cannot find flux file\n\t%s",(*spec)[i].file);
    if (par->thar==1 && access(UVES_replace_envinstr((*spec)[i].thfile),R_OK))
      errormsg("UVES_rinputfile(): Cannot find ThAr file\n\t%s",(*spec)[i].thfile);
    if (par->thar<=1) {
      if (access(UVES_replace_envinstr((*spec)[i].erfile),R_OK))
	errormsg("UVES_rinputfile(): Cannot find error array file\n\t%s\n\
\tfor flux array file\n\t%s",(*spec)[i].erfile,(*spec)[i].file);
    }
    if (access(UVES_replace_envinstr((*spec)[i].wlfile),R_OK))
      errormsg("UVES_rinputfile(): Cannot find wavelen. poly. file\n\t%s\n\
\tfor flux array file\n\t%s",(*spec)[i].wlfile,(*spec)[i].file);
  }

  return 1;

}
