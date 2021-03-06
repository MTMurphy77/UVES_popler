/****************************************************************************
* Read in file specifying observer-frame wavelength regions to mask
* out, usually for atmospheric/telluric feature masking.
****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "UVES_popler.h"
#include "memory.h"
#include "file.h"
#include "error.h"

#define READ_DATA_FILE \
  if ((cptr=fgets(buffer,LNGSTRLEN,data_file))==NULL) { fclose(data_file); \
    errormsg("UVES_ratmask(): Error reading line %d in file %s",i,infile); }
#define ERR_FMT \
  errormsg("UVES_ratmask(): Incorrect format in line %d\n\tof file %s",i,infile);

int UVES_ratmask(atmask *amsk, params *par) {

  int    i=0,j=0,nline=0,ncom=0;
  char   infile[NAMELEN]="\0",buffer[LNGSTRLEN]="\0",dummy[LNGSTRLEN]="\0";
  char   *cptr;
  FILE   *data_file=NULL;

  /* Copy infile name from parameter list structure */
  sprintf(infile,"%s",par->atmaskfile);
  /* Open input file */
  if ((data_file=faskropen("Valid atmospheric mask input file?",infile,5))==NULL)
    errormsg("UVES_ratmask(): Can not open file %s",infile);
  /* Read in first line  */
  i++; READ_DATA_FILE;
  /* Test to see if it's a comment */
  if (sscanf(buffer,"#%s",dummy)==1) ncom++;
  /* Find number of lines in file */
  i++; while ((cptr=fgets(buffer,LNGSTRLEN,data_file))!=NULL) {
    i++; if (sscanf(buffer,"#%s",dummy)==1) ncom++;
  }
  if (!feof(data_file)) {
    fclose(data_file);
    errormsg("UVES_ratmask(): Problem reading line %d in file\n\t%s",i,infile);
  } else { rewind(data_file); nline=i-1; }
  /* Now that the file was successfully read and number of mask
     features is known, transfer properties to atmask structure and
     allocate appropriate memory */
  sprintf(amsk->atmaskfile,"%s",par->atmaskfile);
  amsk->nmask=nline-ncom;
  if ((amsk->swl=darray(amsk->nmask))==NULL)
    errormsg("UVES_ratmask(): Cannot allocate memory for swl\n\
\tarray of size %d for atmospheric mask from file\n\t\%s",amsk->nmask,amsk->atmaskfile);
  if ((amsk->ewl=darray(amsk->nmask))==NULL)
    errormsg("UVES_ratmask(): Cannot allocate memory for ewl\n\
\tarray of size %d for atmospheric mask from file\n\t\%s",amsk->nmask,amsk->atmaskfile);
  /* Read in data  */
  for (i=0,j=0; i<nline; i++) {
    READ_DATA_FILE;
    if (sscanf(buffer,"#%s",dummy)!=1) {
      if (sscanf(buffer,"%lf %lf",&(amsk->swl[j]),&(amsk->ewl[j]))!=2) {
	if (sscanf(buffer,"%lf %lf %s",&(amsk->swl[j]),&(amsk->ewl[j]),dummy)!=3) {
	  fclose(data_file); ERR_FMT;
	}
      }
      j++;
    }
  }
  /* Close file */
  fclose(data_file);

  return 1;

}
