/****************************************************************************
* Read in a 2D FITS echelle spectrum that has been reduced with
* HIRES_REDUX and processed by the KODIAQ team
****************************************************************************/

#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <longnam.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

int UVES_r2Dspec_KODIAQ(spectrum *spec, params *par) {

  double   nulval=0.0;
  double   dwl=0.0;
  long     naxes[9]={0,0,0,0,0,0,0,0,0};
  int      npix=0,col=0;
  int      hdutype=0,hdunum=0,status=0,first=1,naxis=0,anynul=0;
  int      i=0,j=0,k=0;
  int      nphys_ordr=0,phys_ordr_strt=0,*phys_ordr=NULL;
  char     *obj[1],strnul[NAMELEN];
  char     *cptr;
  fitsfile *infits;
  extern  char   *progname;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,UVES_replace_envinstr(spec->file),READONLY,&status))
      errormsg("UVES_r2Dspec_KODIAQ(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  if (fits_get_hdu_type(infits,&hdutype,&status))
    errormsg("UVES_r2Dspec_KODIAQ(): Cannot get HDU type for file\n\t%s",spec->file);
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_KODIAQ(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec_KODIAQ(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (par->thar<=1) {
    if (hdunum<2) errormsg("UVES_r2Dspec_KODIAQ(): Number of HDUs is %d but it\n\
\tmust be either >=%d in file\n\t%s",hdunum,2,spec->file);
  }

  /* Move to second HDU, which contains all the information */
  if (fits_movrel_hdu(infits,1,&hdutype,&status))
    errormsg("UVES_r2Dspec_KODIAQ(): Cannot move to second HDU\n\
\tin FITS file\n\t%s.",spec->file);

  /* Get object name */
  if (par->thar<=1) {
    /* Set the null string value */
    strcpy(strnul," ");
    /* Allocate memory to string array */
    obj[0]=(char *)malloc(NAMELEN);
    /* Get object name column number */
    if (fits_get_colnum(infits,CASEINSEN,"FIELD",&col,&status))
      errormsg("UVES_r2Dspec_KODIAQ(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","FIELD",spec->file);
    /* Read object name */
    if (fits_read_col(infits,TSTRING,col,1,1,1,strnul,&(obj[0]),&anynul,&status)) {
      errormsg("UVES_r2Dspec_KODIAQ(): Cannot read object name\n\
\tfrom column %d of table in FITS file\n\t%s.",1,spec->file);
    }
    sprintf(spec->obj,"%s",obj[0]);
  } else sprintf(spec->obj,"thar_wav");

  /* Alter object name to remove some special characters */
  while ((cptr=strchr(spec->obj,'+'))!=NULL) *cptr='p';
  while ((cptr=strchr(spec->obj,'-'))!=NULL) *cptr='m';
  while ((cptr=strchr(spec->obj,'_'))!=NULL) *cptr='u';
  while ((cptr=strchr(spec->obj,'$'))!=NULL) *cptr='d';
  while ((cptr=strchr(spec->obj,' '))!=NULL) *cptr='s';
  strlower(spec->obj);

  /* Set archival filename */
  sprintf(spec->arfile,"%s",spec->abfile);

  /* Determine the arm we're using, the nominal central wavelength,
     the slit-width, CCD binning, spatial pixel scale, the temperature, 
     atmospheric pressure, observation date */
  /* At the moment, these values are just set to arbitrary numbers */
  spec->cwl=spec->sw=spec->pixscal=spec->temp=spec->pres=-1.0;
  spec->jd=spec->ut=-1.0;
  spec->year=spec->month=spec->day=spec->binx=spec->biny=-1;

  /* Get latitude, longitude and altitude of observatory */
  /* At the moment, these are hard-coded assuming that we really are
     dealing with Keck/HIRES data. Also, they are not really used
     because KODIAQ spectra are already converted to the heliocentric
     frame */
  spec->lat=19.8259; spec->lon=155.4751; spec->alt=4160.0;
  /* Get object RA and DEC and equinox */
  /* These are not available in KODIAQ headers and not really used
     anyway, so just set to arbitrary values here */
  spec->ra=spec->dec=-1.0;

  /* Get exposure time */
  /* Not really needed for processing KODIAQ files so setting this to
     unity */
  spec->etime=1.0;

  /* Read heliocentric velicity correction */
  /* Note: HIRES REDUX files are assumed to have already been
     corrected to the vacuum-heliocentric frame. KODIAQ files are
     normally combinations of multiple exposures, so it is not
     possible to mask out atmospheric lines in general. Therefore,
     this is just set to 0.0 here */
  spec->vhel=0.0;

  /* Read in PHYS_ORDR array to define which orders have valid data */
  if (fits_get_colnum(infits,CASEINSEN,"PHYS_ORDR",&col,&status))
    errormsg("UVES_r2Dspec_KODIAQ(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","PHYS_ORDR",spec->file);
  if (fits_read_tdim(infits,col,9,&naxis,naxes,&status))
    errormsg("UVES_r2Dspec_KODIAQ(): Cannot get dimensions of column named '%s'\n\
\tin binary table in FITS file\n\t%s","PHYS_ORDR",spec->file);
  nphys_ordr=naxis; if (naxis==1) nphys_ordr=naxes[0];
  if ((phys_ordr=iarray(nphys_ordr))==NULL)
    errormsg("UVES_r2Dspec_KODIAQ(): Cannot allocate memory for phys_ordr\n\
\tarray of size %d for file\n\t%s",nphys_ordr,spec->file);
  if (fits_read_col(infits,TINT,col,1,1,nphys_ordr,&nulval,phys_ordr,&anynul,&status))
    errormsg("UVES_r2Dspec_KODIAQ(): Cannot read column %d, '%s', of table\n\
\tin FITS file\n\t%s.",col,"PHYS_ORDR",spec->file);

  /* Determine starting order and number of orders with valid data */
  for (i=0; i<nphys_ordr; i++) {
    if (phys_ordr[i]>0) {
      spec->nor++;
      if (phys_ordr_strt==0) phys_ordr_strt=i;
    }
  }

  /* Allocate memory for echelle order array */
  if (!(spec->or=(echorder *)malloc((size_t)(spec->nor*sizeof(echorder)))))
    errormsg("UVES_r2Dspec_KODIAQ(): Could not allocate memory for echelle\n\
\torder array of size %d",spec->nor);

  /* Find number of pixels to read in for each order */
  if (fits_get_colnum(infits,CASEINSEN,"WAVE",&col,&status))
    errormsg("UVES_r2Dspec_KODIAQ(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","WAVE",spec->file);
  if (fits_read_tdim(infits,col,9,&naxis,naxes,&status))
    errormsg("UVES_r2Dspec_KODIAQ(): Cannot get dimensions of column named '%s'\n\
\tin binary table in FITS file\n\t%s","WAVE",spec->file);
  if (naxis!=2) errormsg("UVES_r2Dspec_KODIAQ(): Expected 2-dimensional column\n\
\tfor '%s' but found a %d-dimensional column instead\n\
\tin binary table in FITS file\n\t%s","WAVE",naxis,spec->file);
  npix=spec->or[0].np=naxes[0];
  for (i=1; i<spec->nor; i++) spec->or[i].np=spec->or[0].np;

  /* Allocate memory for data arrays and fill wavelength, flux, error arrays */
  if (!UVES_memspec(spec,par,spec->nor,0))
    errormsg("UVES_r2Dspec_KODIAQ(): Error returned from UVES_memspec() when\n\
\tattempting to allocate memory for data arrays for file\n\t%s",spec->file);

  /* Read in wavelength information. Note that we read the bluest
     (i.e. highest phys_ordr) order first. Note also that "col" is
     already set to WAVE col number from previous fits_get_colnum
     instance */
  for (i=spec->nor-1,first=phys_ordr_strt*npix+1; i>=0; i--,first+=npix) {
    if (fits_read_col(infits,TDOUBLE,col,1,first,spec->or[i].np,&nulval,
		      spec->or[i].wl,&anynul,&status))
      errormsg("UVES_r2Dspec_KODIAQ(): Cannot read column %d, '%s', of table\n\
\tin FITS file\n\t%s.",col,"WAVE",spec->file);
  }

  /* Read in flux information in same way as wavelength information */
  if (par->thar<=1) {
    if (fits_get_colnum(infits,CASEINSEN,"FX",&col,&status))
      errormsg("UVES_r2Dspec_KODIAQ(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","FX",spec->file);
    for (i=spec->nor-1,first=phys_ordr_strt*npix+1; i>=0; i--,first+=npix) {
      if (fits_read_col(infits,TDOUBLE,col,1,first,spec->or[i].np,&nulval,
			spec->or[i].fl,&anynul,&status))
	errormsg("UVES_r2Dspec_KODIAQ(): Cannot read column %d, '%s', of table\n\
\tin FITS file\n\t%s.",col,"FX",spec->file);
    }
  }

  /* Read in variance information in same way as wavelength
     information and convert it to sigma */
  if (par->thar<=1) {
    if (fits_get_colnum(infits,CASEINSEN,"VAR",&col,&status))
      errormsg("UVES_r2Dspec_KODIAQ(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","VAR",spec->file);
    for (i=spec->nor-1,first=phys_ordr_strt*npix+1; i>=0; i--,first+=npix) {
      if (fits_read_col(infits,TDOUBLE,col,1,first,spec->or[i].np,&nulval,
			spec->or[i].er,&anynul,&status))
	errormsg("UVES_r2Dspec_KODIAQ(): Cannot read column %d, '%s', of table\n\
\tin FITS file\n\t%s.",col,"VAR",spec->file);
      for (j=0; j<spec->or[i].np; j++)
	spec->or[i].er[j]=(spec->or[i].er[j]>0.0) ? sqrt(spec->or[i].er[j]) : -INFIN;
    }
  }

  /* Must treat the special case for KODIAQ files where wavelengths
     are set to zero if there's no information from the chip at those
     wavelengths. At the moment, the fix is just to extrapolate the
     wavelength scale from valid sections of each order out to the
     order edges, only replacing zero wavelength entries. Also,
     sometimes it appears that the dispersion of the wavelength scale
     changes wildly (even sign!) near the order edges, so we chop off
     10 pixels from each edge here. */
  for (i=0; i<spec->nor; i++) {
    /* Find first non-zero wavelength element */
    for (j=0; j<spec->or[i].np-1; j++) if (spec->or[i].wl[j]>DRNDTOL) break;
    dwl=spec->or[i].wl[j+1]-spec->or[i].wl[j];
    for (k=j-1; k>=0; k--) spec->or[i].wl[k]=spec->or[i].wl[k+1]-dwl;
    for (k=j+9; k>=0; k--) spec->or[i].er[k]=-INFIN;
    /* Find last non-zero wavelength element */
    for (j=spec->or[i].np-1; j>=1; j--) if (spec->or[i].wl[j]>DRNDTOL) break;
    dwl=spec->or[i].wl[j]-spec->or[i].wl[j-1];
    for (k=j+1; k<spec->or[i].np; k++) spec->or[i].wl[k]=spec->or[i].wl[k-1]+dwl;
    for (k=j-9; k<spec->or[i].np; k++) spec->or[i].er[k]=-INFIN;
  }

  /* Close flux fits file */
  fits_close_file(infits,&status);

  /* If ThAr information is to be read in, do it now */
  if (par->thar==1) {
    /* Actually, this is not implemented for KODIAQ files. Sorry. */
     errormsg("UVES_r2Dspec_KODIAQ(): Reading KODIAQ FITS is not currently\n\
\tpossible in %s. Sorry!",progname);
  }

  /* Set order identity numbers and parameters */
  for (i=0; i<spec->nor; i++) {
    spec->or[i].id=i; spec->or[i].sid=spec->id;
    /* Find the number of useful pixels in the order */
    if (par->thar<=1) {
      for (spec->or[i].nuse=0,j=0; j<spec->or[i].np; j++) {
	if (spec->or[i].er[j]>DRNDTOL) spec->or[i].nuse++;
	else spec->or[i].st[j]=RCLIP;
	if (par->thar==1 && spec->or[i].th[j]==0.0) spec->or[i].tst[j]=RCLIP;
      }
    } else {
      /* Mask bad pixels */
      for (j=0; j<spec->or[i].np; j++)
	if (spec->or[i].fl[j]<=0.0) spec->or[i].st[j]=RCLIP;
      /* Determine where the useful part of the spectrum begins and ends */
      j=0; while (j<spec->or[i].np && spec->or[i].st[j]<0) j++;
      k=j; j=spec->or[i].np-1; while (j>=0 && spec->or[i].st[j]<0) j--;
      spec->or[i].nuse=j-k+1;
    }
    if (spec->or[i].nuse<MINUSE) {
      free(spec->or[i].wl);
      free(spec->or[i].vhwl); free(spec->or[i].vhrwl); free(spec->or[i].fl);
      free(spec->or[i].er); free(spec->or[i].res); free(spec->or[i].st);
      if (par->thar==1) {
	free(spec->or[i].th); free(spec->or[i].ter); free(spec->or[i].tst);
      }
    }
  }

  /* Set median resolution to zero for lack of any other information */
  spec->arcfwhm=0.0;

  /* Make sure this spectrum will be combined into combined spectrum later */
  spec->comb=1;

  return 1;

}
