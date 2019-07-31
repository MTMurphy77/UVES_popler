/****************************************************************************
* Read in a UVES 2D FITS echelle spectrum and error array
****************************************************************************/

#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <longnam.h>
#include "UVES_popler.h"
#include "stats.h"
#include "memory.h"
#include "error.h"
#include "const.h"

#define FITS_RKEYS(KEYWORD) { \
  if (strlen(KEYWORD)>1) sprintf(key,"%s",KEYWORD); sprintf(cards,"%s","\0"); \
  if (fits_read_key(infits,TSTRING,key,cards,comment,&status)) \
    errormsg("UVES_r2Dspec_espresso(): Cannot read header keyword\n\t%s in FITS file\n\t%s.", \
	     key,spec->file); }
#define FITS_RKEYD(KEYWORD) { \
  if (strlen(KEYWORD)>1) sprintf(key,"%s",KEYWORD); \
  if (fits_read_key(infits,TDOUBLE,key,&cardd,comment,&status)) \
    errormsg("UVES_r2Dspec_espresso(): Cannot read header keyword\n\t%s in FITS file\n\t%s.", \
	     key,spec->file); }
#define FITS_RKEYI(KEYWORD) { \
  if (strlen(KEYWORD)>1) sprintf(key,"%s",KEYWORD); \
  if (fits_read_key(infits,TINT,key,&cardi,comment,&status)) \
    errormsg("UVES_r2Dspec_espresso(): Cannot read header keyword\n\t%s in FITS file\n\t%s.", \
	     key,spec->file); }

int UVES_r2Dspec_espresso(spectrum *spec, params *par) {

  double   cardd=0.0;
  double   hour=0.0,min=0.0,sec=0.0;
  double   nulval=0.0;
  long     naxes[9]={0,0,0,0,0,0,0,0,0};
  int      cardi=0;
  int      npix=0,telnum=0;
  int      hdutype=0,hdunum=0,status=0,bitpix=0,first=1,naxis=0,anynul=0;
  int      i=0,j=0,k=0;
  char     comment[FLEN_COMMENT]="\0";
  char     key[FLEN_KEYWORD]="\0";
  char     date[FLEN_CARD]="\0",instrument[FLEN_CARD]="\0";
  char     cards[FLEN_CARD]="\0";
  char     *cptr;
  fitsfile *infits;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,UVES_replace_envinstr(spec->file),READONLY,&status))
      errormsg("UVES_r2Dspec_espresso(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  if (fits_get_hdu_type(infits,&hdutype,&status))
    errormsg("UVES_r2Dspec_espresso(): Cannot get HDU type for file\n\t%s",spec->file);
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_espresso(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec_espresso(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (hdunum<5) errormsg("UVES_r2Dspec_espresso(): Number of HDUs is %d instead\n\
\tof %d (at least) in file\n\t%s",hdunum,5,spec->file);

  /* Check that this is really an ESPRESSO file */
  FITS_RKEYS("INSTRUME"); sprintf(instrument,"%s",cards);
  if (!strncmp(instrument,"ESPRESSO",8)) spec->fvers=0;
  else errormsg("UVES_r2Dspec_espresso(): Header keyword %s='%s' is not\n\
\t'ESPRESSO'. This file,\n\t%s\n\
\tdoesn't look like it's from ESPRESSO.",key,cards,spec->file);

  /* Read in entire main header as array of strings */
  if (fits_get_hdrspace(infits,&(spec->nhead_ori),NULL,&status))
      errormsg("UVES_r2Dspec_espresso(): Cannot read number of keys in main\n\
\theader in FITS file\n\t%s.",spec->file);
  if ((spec->head_ori=cmatrix(spec->nhead_ori,FLEN_CARD))==NULL)
    errormsg("UVES_r2Dspec_espresso(): Cannot allocate memory for head_ori\n\
\tmatrix of size %dx%d for file %s\n",spec->nhead_ori,FLEN_CARD,spec->file);
  for (i=0; i<spec->nhead_ori; i++) {
    if (fits_read_record(infits,i+1,spec->head_ori[i],&status))
      errormsg("UVES_r2Dspec_espresso(): Cannot read header key %d from\n\
\tmain header of file %s\n",i+1,spec->file);
  }

  /* Get object name */
  if (par->thar<=1) {
    sprintf(key,"HIERARCH ESO OBS TARG NAME"); FITS_RKEYS("");
    sprintf(spec->obj,"%s",cards);
  } else sprintf(spec->obj,"thar_wav");
  /* Alter object name to remove some special characters and convert
     to lower case*/
  while ((cptr=strchr(spec->obj,'+'))!=NULL) *cptr='p';
  while ((cptr=strchr(spec->obj,'-'))!=NULL) *cptr='m';
  while ((cptr=strchr(spec->obj,'_'))!=NULL) *cptr='u';
  while ((cptr=strchr(spec->obj,'$'))!=NULL) *cptr='d';
  while ((cptr=strchr(spec->obj,' '))!=NULL) *cptr='s';
  strlower(spec->obj);

  /* Read in archival filename */
  if (!spec->fvers) sprintf(spec->arfile,"%s",spec->file);
  else { FITS_RKEYS("FILENAME"); sprintf(spec->arfile,"%s",cards); }

  /* Read in type of observation, i.e. whether we should expect
     another exposure containing the sky */
  /* FOR THE MOMENT, WE ASSUME REDUCED ESPRESSO PRODUCTS ARE EITHER
     SKY-SUBTRACT, AS SUPPLIED, OR THAT NO SKY FIBRE WAS USED */
  spec->skysub=0;

  /* Determine the arm we're using, the nominal central wavelength,
     the slit-width, CCD binning, spatial pixel scale, the temperature
     and the atmospheric pressure */
  /* At the moment, some of these values are just set to arbitrary numbers */
  spec->cwl=spec->sw=spec->pixscal=-1.0;
  spec->binx=spec->biny=-1;
  /* This should be the temperature at the grating */
  if (!spec->fvers) {
    sprintf(key,"HIERARCH ESO INS5 TEMP10 VAL"); FITS_RKEYD("");
    spec->temp=cardd;
    /* Note that ESPRESSO is a vacuum spectrograph so the atmospheric
       pressure read in below is not relevant to the spectrograph itself */
    telnum=1; while (telnum<4) {
      sprintf(key,"HIERARCH ESO TEL%d AMBI PRES START",telnum);
      if (fits_read_key(infits,TDOUBLE,key,&cardd,comment,&status)) { telnum++; status=0; }
      else break;
    }
    if (telnum==4) errormsg("UVES_r2Dspec_espresso(): Cannot read header keyword\n\
\t%s, for ?=1-4, in FITS file\n\t%s.","HIERARCH ESO TEL? AMBI PRES START",spec->file);
    spec->pres=cardd;
    sprintf(key,"HIERARCH ESO TEL%d AMBI PRES END",telnum); FITS_RKEYD("");
    spec->pres+=cardd; spec->pres*=0.5;
    FITS_RKEYD("HIERARCH ESO QC BERV"); spec->vbar_head=cardd;
  }

  if (par->thar<=1) {
    /* Get modified julian day */
    FITS_RKEYD("MJD-OBS"); spec->jd=cardd;
    /* Convert to Julian day */
    spec->jd+=2400000.5;
    /* Get date of observation and convert to year, month and day */
    FITS_RKEYS("DATE-OBS"); sprintf(date,"%s",cards);
    if (strlen(date)>10) {
      /* This means that the UT is appended to the date */
      if (sscanf(date,"%d-%d-%dT%lf:%lf:%lf",&(spec->year),&(spec->month),
		 &(spec->day),&hour,&min,&sec)!=6)
	errormsg("UVES_r2Dspec_espresso(): Cannot read format of keyword\n\
%s='%s' in FITS file\n\t%s",key,date,spec->file);
      /* Convert hours, mins and secs into the UT */
      spec->ut=sec/3600.0+min/60.0+hour;
    } else errormsg("UVES_r2Dspec_espresso(): Cannot read format of keyword\n\
%s='%s' in FITS file\n\t%s",key,date,spec->file);

    /* Get latitude, longitude and altitude of observatory */
    if (!spec->fvers) {
      telnum=1; while (telnum<4){
	sprintf(key,"HIERARCH ESO TEL%d GEOLAT",telnum);
	if (fits_read_key(infits,TDOUBLE,key,&cardd,comment,&status)) { telnum++; status=0; }
	else break;
      }
      if (telnum==4) errormsg("UVES_r2Dspec_espresso(): Cannot read header keyword\n\
\t%s, for ?=1-4, in FITS file\n\t%s.","HIERARCH ESO TEL? GEOLAT",spec->file);
      spec->lat=cardd;
      sprintf(key,"HIERARCH ESO TEL%d GEOLON",telnum); FITS_RKEYD(""); spec->lon=cardd;
      /* Convert Longitude to normal convention. ESO files use East as
	 positive longitudes. However, the normal convention is for West
	 to be positive. */
      if (spec->lon<0.0) spec->lon*=-1.0;
      else spec->lon=360.0-spec->lon;
      sprintf(key,"HIERARCH ESO TEL%d GEOELEV",telnum); FITS_RKEYD("");  spec->alt=cardd;
    }

    /* Get object RA and DEC and equinox */
    if (!spec->fvers) { FITS_RKEYD("RA"); }
    spec->ra=cardd;
    /* Convert RA into hours */
    spec->ra/=15.0;
    if (!spec->fvers) { FITS_RKEYD("DEC"); }
    spec->dec=cardd;
    FITS_RKEYD("EQUINOX"); spec->equ=cardd;
  }

  /* Get exposure time */
  FITS_RKEYD("EXPTIME"); spec->etime=cardd;


  /* Move to second HDU to read in flux */
  if (!spec->fvers && fits_movrel_hdu(infits,1,&hdutype,&status))
    errormsg("UVES_r2Dspec_espresso(): Could not move to second HDU\n\
\tin file %s",spec->file);

  /* Check that this is an IMAGE HDU */
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_espresso(): Second HDU is not a FITS image: %s",spec->file);

  /* Get number of echelle orders & allocate memory for echelle order array */
  FITS_RKEYI("NAXIS2"); spec->nor=cardi;
  if (!(spec->or=(echorder *)malloc((size_t)(spec->nor*sizeof(echorder)))))
    errormsg("UVES_r2Dspec_espresso(): Could not allocate memory for echelle\n\
\torder array of size %d for spectrum in file\n\t%s",spec->nor,spec->file);

  /* Find number of pixels to read in for each order */
  FITS_RKEYI("NAXIS1"); spec->or[0].np=cardi;
  for (i=1; i<spec->nor; i++) spec->or[i].np=spec->or[0].np;

  /* Allocate memory for data arrays */
  if (!UVES_memspec(spec,par,spec->nor,0))
    errormsg("UVES_r2Dspec_espresso(): Error returned from UVES_memspec() when\n\
\tattempting to allocate memory for data arrays for file\n\t%s",spec->file);

  /* Get image dimensions */
  if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
    errormsg("UVES_r2Dspec_espresso(): Couldn't get image dimensions for\n\
\tfile %s",spec->file);
  if (!naxis) errormsg("UVES_r2Dspec_espresso(): Couldn't find input image extension\n\
\tfor file %s",spec->file);
  npix=naxes[0];

  /* Read in flux information */
  for (i=0,first=1; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].fl,&anynul,&status))
      errormsg("UVES_r2Dspec_espresso(): Cannot read flux array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }


  /* Move to third HDU to read in error */
  if (!spec->fvers && fits_movrel_hdu(infits,1,&hdutype,&status))
    errormsg("UVES_r2Dspec_espresso(): Could not move to third HDU\n\
\tin file %s",spec->file);

  /* Check that this is an IMAGE HDU */
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_espresso(): Third HDU is not a FITS image: %s",spec->file);

  /* Check number of echelle orders */
  FITS_RKEYI("NAXIS2");
  if (cardi!=spec->nor) errormsg("UVES_r2Dspec_espresso(): Number of orders in error HDU (%d)\n\
\tis not equal to number in flux HDU (%d) in file\n\t%s",cardi,spec->nor,spec->file);

  /* Find number of pixels to read in for each order */
  FITS_RKEYI("NAXIS1");
  if (cardi!=spec->or[0].np) errormsg("UVES_r2Dspec_espresso(): Number of pixels per order\n\
\tin error HDU (%d) is not equal to number in flux HDU (%d) in file\n\t%s",
				 cardi,spec->or[0].np,spec->file);

  /* Read in error information */
  for (i=0,first=1; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].er,&anynul,&status))
      errormsg("UVES_r2Dspec_espresso(): Cannot read error array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }


  /* Move to third HDU to read in status */
  if (!spec->fvers && fits_movrel_hdu(infits,1,&hdutype,&status))
    errormsg("UVES_r2Dspec_espresso(): Could not move to third HDU\n\
\tin file %s",spec->file);

  /* Check that this is an IMAGE HDU */
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_espresso(): Third HDU is not a FITS image: %s",spec->file);

  /* Check number of echelle orders */
  FITS_RKEYI("NAXIS2");
  if (cardi!=spec->nor) errormsg("UVES_r2Dspec_espresso(): Number of orders in status HDU (%d)\n\
\tis not equal to number in flux HDU (%d) in file\n\t%s",cardi,spec->nor,spec->file);

  /* Find number of pixels to read in for each order */
  FITS_RKEYI("NAXIS1");
  if (cardi!=spec->or[0].np) errormsg("UVES_r2Dspec_espresso(): Number of pixels per order\n\
\tin status HDU (%d) is not equal to number in flux HDU (%d) in file\n\t%s",
				 cardi,spec->or[0].np,spec->file);

  /* Read in status information */
  for (i=0,first=1; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].wl,&anynul,&status))
      errormsg("UVES_r2Dspec_espresso(): Cannot read status array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
    for (j=0; j<spec->or[i].np; j++) spec->or[i].st[j]=(spec->or[i].wl[j]>0.0) ? RCLIP : 1;
  }


  /* Move to fourth HDU to read in wavelength */
  if (!spec->fvers && fits_movrel_hdu(infits,1,&hdutype,&status))
    errormsg("UVES_r2Dspec_espresso(): Could not move to fourth HDU\n\
\tin file %s",spec->file);

  /* Check that this is an IMAGE HDU */
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_espresso(): Fourth HDU is not a FITS image: %s",spec->file);

  /* Check number of echelle orders */
  FITS_RKEYI("NAXIS2");
  if (cardi!=spec->nor) errormsg("UVES_r2Dspec_espresso(): Number of orders in wavelength HDU (%d)\n\
\tis not equal to number in flux HDU (%d) in file\n\t%s",cardi,spec->nor,spec->file);

  /* Find number of pixels to read in for each order */
  FITS_RKEYI("NAXIS1");
  if (cardi!=spec->or[0].np) errormsg("UVES_r2Dspec_espresso(): Number of pixels per order\n\
\tin wavelength HDU (%d) is not equal to number in flux HDU (%d) in file\n\t%s",
				 cardi,spec->or[0].np,spec->file);

  /* Read in wavelength information */
  for (i=0,first=1; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].wl,&anynul,&status))
      errormsg("UVES_r2Dspec_espresso(): Cannot read wavelength array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }

  /* Close flux fits file */
  fits_close_file(infits,&status);


  /* If ThAr information is to be read in, do it now */
  if (par->thar==1) {
    errormsg("UVES_r2Dspec_espresso(): Currently there is no functionality\n\
\tto read in the ThAr spectrum from ESPRESSO files. Sorry.");
  }


  /* Set order identity numbers and parameters in flux spectrum */
  for (i=0; i<spec->nor; i++) {
    spec->or[i].id=i; spec->or[i].sid=spec->id;
    /* Find the number of useful pixels in the order */
    if (par->thar<=1) {
      for (spec->or[i].nuse=0,j=0; j<spec->or[i].np; j++) {
	if (spec->or[i].er[j]>0.0 && spec->or[i].st[j]==1) spec->or[i].nuse++;
	else spec->or[i].st[j]=RCLIP;
	if (par->thar==1) {
	  if (spec->or[i].th[j]==0.0) spec->or[i].tst[j]=RCLIP;
	}
      }
    } else {
      /* Mask bad pixels */
      for (j=0; j<spec->or[i].np; j++)
	if (spec->or[i].fl[j]==0.0) spec->or[i].st[j]=RCLIP;
      /* Determine where the useful part of the spectrum begins and ends */
      j=0; while (j<spec->or[i].np && spec->or[i].st[j]<0) j++;
      k=j; j=spec->or[i].np-1; while (j>=0 && spec->or[i].st[j]<0) j--;
      spec->or[i].nuse=j-k+1;
    }
    if (spec->or[i].nuse<MINUSE) {
      /* Free the allocated memory for the raw arrays */
      if (!UVES_memspec(spec,par,-i-1,0))
	errormsg("UVES_r2Dspec_espresso(): Error returned from UVES_memspec() when\n\
\tattempting to free memory from order %d of file\n\t%s",i+1,spec->file);
    }
  }

  /* Make sure this spectrum will be combined into combined spectrum later */
  spec->comb=1;

  return 1;

}
