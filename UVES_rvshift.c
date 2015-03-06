/****************************************************************************
* Read in file specifying velocity shifts (plus higher-order
* distortions) to be applied to individual spectra.
****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "UVES_popler.h"
#include "file.h"
#include "error.h"

#define READ_DATA_FILE \
  if ((cptr=fgets(buffer,LNGSTRLEN,data_file))==NULL) { fclose(data_file); \
    errormsg("UVES_rvshift(): Error reading line %d in file %s",i,infile); }
#define ERR_FMT \
  errormsg("UVES_rvshift(): Incorrect format in line %d\n\tof file %s",i+1,infile);

int UVES_rvshift(spectrum **spec, int nspec, params *par) {

  double vshift=0.0,vslope=0.0,refwav=0.0;
  int    i=0,j=0,nline=0,ninput=0;
  char   infile[NAMELEN]="\0",buffer[LNGSTRLEN]="\0";
  char   *cptr;
  char   vshiftfile[LNGSTRLEN];
  FILE   *data_file=NULL;

  /* Copy infile name from parameter list structure */
  sprintf(infile,"%s",par->vshiftfile);
  /* Open input file */
  if ((data_file=faskropen("Valid velocity shift input file?",infile,5))==NULL)
    errormsg("UVES_rvshift(): Can not open file %s",infile);
  /* Read in first line  */
  i++; READ_DATA_FILE;
  /* Find number of lines in file */
  i++; while ((cptr=fgets(buffer,LNGSTRLEN,data_file))!=NULL) i++;
  if (!feof(data_file)) {
    fclose(data_file);
    errormsg("UVES_rvshift(): Problem reading line %d in file\n\t%s",i,infile);
  } else { rewind(data_file); nline=i-1; }
  /* Read in data and try to match file names with those of the actual
     spectra. If they match, record velocity shift and distortion
     information in each spectrum's structure */
  for (i=0; i<nline; i++) {
    READ_DATA_FILE;
    vshift=vslope=refwav=0.0;
    if ((ninput=sscanf(buffer,"%s %lf %lf %lf",vshiftfile,&vshift,&vslope,&refwav))!=4) {
      if (ninput!=2) { fclose(data_file); ERR_FMT; }
      else sscanf(buffer,"%s %lf",vshiftfile,&vshift);
    }
    for (j=0; j<nspec; j++) if (!strcmp(vshiftfile,(*spec)[j].file)) break;
    if (j==nspec)
      errormsg("UVES_rvshift(): Cannot find file\n\t%s\n\
\tspecified in velocity shift file\n\t%s\n\tin list of flux files");
    (*spec)[j].vshift=vshift; (*spec)[j].vslope=vslope; (*spec)[j].refwav=refwav;
  }
  /* Close file */
  fclose(data_file);

  return 1;

}
