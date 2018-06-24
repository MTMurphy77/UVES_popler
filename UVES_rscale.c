/****************************************************************************
* Read in file specifying flux/error scalings to be applied to
* individual spectra, orders and/or the combined spectrum
****************************************************************************/

#include <string.h>
#include "UVES_popler.h"
#include "file.h"
#include "error.h"

#define READ_DATA_FILE \
  if ((cptr=fgets(buffer,LNGSTRLEN,data_file))==NULL) { fclose(data_file); \
    errormsg("UVES_rscale(): Error reading line %d in file %s",i,infile); }
#define ERR_FMT \
  errormsg("UVES_rscale(): Incorrect format in line %d\n\tof file %s",i+1,infile);

int UVES_rscale(spectrum **spec, int nspec, cspectrum *cspec, params *par) {

  double scale=0.0;
  int    i=0,j=0,k=0,sord=0,eord=0,insclfe=-1,insclbefaft=0,nline=0,ninput=0;
  char   infile[NAMELEN]="\0",buffer[LNGSTRLEN]="\0";
  char   *cptr;
  char   scalefile[LNGSTRLEN];
  char   ordstr[NAMELEN]="\0",flxerrstr[NAMELEN]="\0",befaftstr[NAMELEN]="\0";
  char   insclstr[VLNGSTRLEN]="\0";
  FILE   *data_file=NULL;

  /* Copy infile name from parameter list structure */
  sprintf(infile,"%s",par->scalefile);
  /* Open input file */
  if ((data_file=faskropen("Valid flux/error scaling input file?",infile,5))==NULL)
    errormsg("UVES_rscale(): Can not open file %s",infile);
  /* Read in first line  */
  i++; READ_DATA_FILE;
  /* Find number of lines in file */
  i++; while ((cptr=fgets(buffer,LNGSTRLEN,data_file))!=NULL) i++;
  if (!feof(data_file)) {
    fclose(data_file);
    errormsg("UVES_rscale(): Problem reading line %d in file\n\t%s",i,infile);
  } else { rewind(data_file); nline=i-1; }
  /* Read in data and try to match file names with those of the actual
     spectra. If they match, record scaling information in the
     relevant spectrum's structure or the combined spectrum's
     structure */
  for (i=0; i<nline; i++) {
    READ_DATA_FILE;
    if ((ninput=sscanf(buffer,"%s o%s %s %lf %s",scalefile,ordstr,flxerrstr,&scale,
		       befaftstr))!=5) { fclose(data_file); ERR_FMT; }
    for (j=0; j<nspec; j++) if (strstr((*spec)[j].file,scalefile)) break;
    if (j==nspec) errormsg("UVES_rscale(): Cannot find file\n\t%s\n\
\tspecified in flux/error scaling file\n\t%s\n\
\tin list of flux files",scalefile,infile);
    /* Read which order(s) should be scaled */
    if (!strncmp(ordstr,"all",3)) { sord=-1; eord=-1; }
    else if (sscanf(ordstr,"%d-%d",&sord,&eord)==2) {
      if (sord<0 || sord>eord)
	errormsg("UVES_rscale(): Starting and ending order numbers (=%d-%d)\n\
\tspecified in flux/error scaling list of flux files\n\
\t%s\n\tare not sensible for scaling spectrum in file\n\t%s",
		 sord,eord,infile,scalefile);
    } else if (sscanf(ordstr,"%d",&sord)==1) {
      if (sord>=0) eord=sord;
      else errormsg("UVES_rscale(): Order number (=%d) specified in flux/error\n\
\tscaling list of flux files\n\t%s\n\
\tis not sensible for scaling spectrum in file\n\t%s",sord,infile,scalefile);
    } else errormsg("UVES_rscale(): Could not read format for which orders to\n\
\tscale in file %s\n\tas specified in flux/error scaling file\n\t%s.\n\
\tFormat given is o%s but allowed formats are o#, o#-#, or oall",
		    scalefile,infile,ordstr);
    /* Read whether scaling is for flux, error or both */
    if (strstr(flxerrstr,"f")) insclfe=0;
    if (strstr(flxerrstr,"e")) insclfe=(insclfe==0) ? 2 : 1;
    if (insclfe==-1)
      errormsg("UVES_rscale(): Could not read format for scaling either flux\n\
\tor error array of orders %d-%d in file\n\t%s\n\
\tas specified in flux/error scaling file\n\t%s.\n\
\tFormat given is %s but allowed formats are f, e, fe or ef",
	       sord,eord,scalefile,infile,flxerrstr);
    /* Read whether scaling is before or after automatic section and manual actions */
    if (strstr(befaftstr,"b")) insclbefaft=1;
    else if (strstr(befaftstr,"a")) insclbefaft=2;
    else errormsg("UVES_rscale(): Could not read format for scaling before or after\n\
\tautomatic & manual actions for orders %d-%d in file\n\t%s\n\
\tas specified in flux/error scaling file\n\t%s.\n\
\tFormat given is %s but allowed formats are b or a",
	       sord,eord,scalefile,infile,befaftstr);
    /* Loop through scaled orders and add to scaling string for this spectrum */
    for (k=sord; k<=eord; k++) {
      sprintf(insclstr,"%d,%d,%lf,%d;",k,insclfe,scale,insclbefaft);
      sprintf((*spec)[j].inscl+strlen((*spec)[j].inscl),"%s",insclstr);
    }
  }
  /* Close file */
  fclose(data_file);

  return 1;

}
