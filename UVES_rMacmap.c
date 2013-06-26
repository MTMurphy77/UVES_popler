/****************************************************************************
* Read in file for mapping path and file names in UPL file to
* user-defined names. This is typically done using a Macmap file,
* written out by UVES_headsort, for mapping case-sensitive
* path/file-names to case-insensitive ones.
****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "UVES_popler.h"
#include "memory.h"
#include "file.h"
#include "error.h"

#define READ_DATA_FILE \
  if ((cptr=fgets(buffer,LNGSTRLEN,data_file))==NULL) { fclose(data_file); \
    errormsg("UVES_rMacmap(): Error reading line %d in file %s",i,infile); }
#define ERR_FMT \
  errormsg("UVES_rMacmap(): Incorrect format in line %d\n\tof file %s",i,infile);

int UVES_rMacmap(macmap *mmap) {

  int    i=0;
  char   infile[NAMELEN]="\0",buffer[LNGSTRLEN]="\0";
  char   *cptr;
  FILE   *data_file=NULL;

  /* Copy infile name from macmap structure */
  sprintf(infile,"%s",mmap->mmapfile);
  /* Open input file */
  if ((data_file=faskropen("Valid Macmap input file?",infile,5))==NULL)
    errormsg("UVES_rMacmap(): Can not open file %s",infile);
  /* Read in first line  */
  i++; READ_DATA_FILE;
  /* Find number of lines in file */
  i++; while ((cptr=fgets(buffer,LNGSTRLEN,data_file))!=NULL) i++;
  if (!feof(data_file)) {
    fclose(data_file);
    errormsg("UVES_rMacmap(): Problem reading line %d in file\n\t%s",i,
	     infile);
  } else { rewind(data_file); mmap->nmap=i-1; }
  /* Allocate memory for Macmap entires */
  if ((mmap->idx_ori=iarray(mmap->nmap))==NULL)
    errormsg("UVES_rMacmap(): Cannot allocate memory for idx_ori\n\
\tarray of length %d for file %s\n",mmap->nmap,infile);
  if ((mmap->idx_new=iarray(mmap->nmap))==NULL)
    errormsg("UVES_rMacmap(): Cannot allocate memory for idx_new\n\
\tarray of length %d for file %s\n",mmap->nmap,infile);
  if ((mmap->obj_ori=cmatrix(mmap->nmap,NAMELEN))==NULL)
    errormsg("UVES_rMacmap(): Cannot allocate memory for obj_ori\n\
\tmatrix of size %dx%d for file %s\n",mmap->nmap,NAMELEN,infile);
  if ((mmap->obj_new=cmatrix(mmap->nmap,NAMELEN))==NULL)
    errormsg("UVES_rMacmap(): Cannot allocate memory for obj_new\n\
\tmatrix of size %dx%d for file %s\n",mmap->nmap,NAMELEN,infile);
  if ((mmap->cwl_ori=cmatrix(mmap->nmap,NAMELEN))==NULL)
    errormsg("UVES_rMacmap(): Cannot allocate memory for cwl_ori\n\
\tmatrix of size %dx%d for file %s\n",mmap->nmap,NAMELEN,infile);
  if ((mmap->cwl_new=cmatrix(mmap->nmap,NAMELEN))==NULL)
    errormsg("UVES_rMacmap(): Cannot allocate memory for cwl_new\n\
\tmatrix of size %dx%d for file %s\n",mmap->nmap,NAMELEN,infile);
  /* Read in data */
  for (i=0; i<mmap->nmap; i++) {
    READ_DATA_FILE;
    if (sscanf(buffer,"%s %s %d %s %s %d",mmap->obj_ori[i],mmap->cwl_ori[i],
	       &(mmap->idx_ori[i]),mmap->obj_new[i],mmap->cwl_new[i],
	       &(mmap->idx_new[i]))!=6) {
      fclose(data_file); ERR_FMT;
    }
  }
  /* Close file */
  fclose(data_file);

  return 1;

}
