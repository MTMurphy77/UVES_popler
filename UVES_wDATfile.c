/****************************************************************************
* Write out a data file whic contains the normalised spectrum and error
****************************************************************************/

#include <stdlib.h>
#include <unistd.h>
#include "UVES_popler.h"
#include "file.h"
#include "error.h"

int UVES_wDATfile(char *filename, cspectrum *cspec) {

  int       i=0;
  FILE      *data_file=NULL;

  /* Open specified file */
  if (!access(filename,R_OK))
    warnmsg("Overwriting existing file %s",filename);
  else fprintf(stderr,"Opening output data file %s\n",filename);
  if ((data_file=faskwopen("Output spectrum data file name?",filename,4))==NULL)
    errormsg("Can not open file %s",filename);

  /* Write out normalized spectrum to file */
  for (i=0; i<cspec->np; i++)
    fprintf(data_file,"%16.10lf %19.14lf %19.14lg %19.14lg %19.14lg\n",
	    cspec->wl[i],cspec->no[i],cspec->ne[i],cspec->nf[i],cspec->co[i]);

  /* Close file */
  fclose(data_file);

  return 1;

}
