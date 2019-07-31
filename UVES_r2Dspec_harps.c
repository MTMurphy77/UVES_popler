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
    errormsg("UVES_r2Dspec_harps(): Cannot read header keyword\n\t%s in FITS file\n\t%s.", \
	     key,spec->file); }
#define FITS_RKEYD(KEYWORD) { \
  if (strlen(KEYWORD)>1) sprintf(key,"%s",KEYWORD); \
  if (fits_read_key(infits,TDOUBLE,key,&cardd,comment,&status)) \
    errormsg("UVES_r2Dspec_harps(): Cannot read header keyword\n\t%s in FITS file\n\t%s.", \
	     key,spec->file); }
#define FITS_RKEYI(KEYWORD) { \
  if (strlen(KEYWORD)>1) sprintf(key,"%s",KEYWORD); \
  if (fits_read_key(infits,TINT,key,&cardi,comment,&status)) \
    errormsg("UVES_r2Dspec_harps(): Cannot read header keyword\n\t%s in FITS file\n\t%s.", \
	     key,spec->file); }

int UVES_r2Dspec_harps(spectrum *spec, params *par) {

  double   cardd=0.0;
  double   hour=0.0,deg=0.0,min=0.0,sec=0.0,eperADU=0.0,rdnoise=0.0;
  double   nulval=0.0;
  double   **blz=NULL;
  long     naxes[9]={0,0,0,0,0,0,0,0,0};
  int      cardi=0,neffspatpix=0;
  int      ndeg=0,npixobj=0,npixsky=0;
  int      hdutype=0,hdunum=0,status=0,bitpix=0,first=1,naxis=0,anynul=0;
  int      i=0,j=0,k=0;
  char     comment[FLEN_COMMENT]="\0";
  char     key[FLEN_KEYWORD]="\0";
  char     date[FLEN_CARD]="\0",observatory[FLEN_CARD]="\0",instrument[FLEN_CARD]="\0";
  char     cards[FLEN_CARD]="\0",hem[FLEN_CARD]="\0";
  char     blzfile[LNGSTRLEN]="\0",sblzfile[LNGSTRLEN]="\0";
  char     *cptr;
  fitsfile *infits;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,UVES_replace_envinstr(spec->file),READONLY,&status))
      errormsg("UVES_r2Dspec_harps(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  if (fits_get_hdu_type(infits,&hdutype,&status))
    errormsg("UVES_r2Dspec_harps(): Cannot get HDU type for file\n\t%s",spec->file);
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_harps(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec_harps(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (hdunum>1) errormsg("UVES_r2Dspec_harps(): Number of HDUs is %d instead\n\
\tof %d (at most) in file\n\t%s",hdunum,1,spec->file);

  /* Check which HARPS we are dealing with here, ESO ("HARPS") or TNG ("HARPN") */
  FITS_RKEYS("INSTRUME"); sprintf(instrument,"%s",cards);
  if (!strncmp(instrument,"HARPS",5)) spec->fvers=0;
  else if (!strncmp(instrument,"HARPN",5)) spec->fvers=1;
  else errormsg("UVES_r2Dspec_harps(): Header keyword %s='%s' is neither\n\
\t'HARPS' nor 'HARPN'. This file,\n\t%s\n\
\tdoesn't look like it's from one of the HARPS instruments.",key,cards,spec->file);
  if (!spec->fvers) sprintf(observatory,"ESO");
  else sprintf(observatory,"TNG");
  /* Get data reduction software version */
  /* sprintf(key,"HIERARCH %s DRS VERSION",observatory); FITS_RKEYS(""); */

  /* Read in entire main header as array of strings */
  if (fits_get_hdrspace(infits,&(spec->nhead_ori),NULL,&status))
      errormsg("UVES_r2Dspec_harps(): Cannot read number of keys in main\n\
\theader in FITS file\n\t%s.",spec->file);
  if ((spec->head_ori=cmatrix(spec->nhead_ori,FLEN_CARD))==NULL)
    errormsg("UVES_r2Dspec_harps(): Cannot allocate memory for head_ori\n\
\tmatrix of size %dx%d for file %s\n",spec->nhead_ori,FLEN_CARD,spec->file);
  for (i=0; i<spec->nhead_ori; i++) {
    if (fits_read_record(infits,i+1,spec->head_ori[i],&status))
      errormsg("UVES_r2Dspec_harps(): Cannot read header key %d from\n\
\tmain header of file %s\n",i+1,spec->file);
  }

  /* Get object name */
  if (par->thar<=1) {
    sprintf(key,"HIERARCH %s OBS TARG NAME",observatory); FITS_RKEYS("");
    sprintf(spec->obj,"%s",cards);
  }
  else sprintf(spec->obj,"thar_wav");
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
  spec->skysub=0;
  sprintf(key,"HIERARCH %s DPR TYPE",observatory); FITS_RKEYS("");
  if (strstr(cards,"SKY")!=NULL) {
    spec->skysub=1; sprintf(spec->ss.path,"%s",spec->path);
    sprintf(spec->ss.abfile,"%s",spec->abfile);
    if ((cptr=strstr(spec->ss.abfile,"e2ds_A"))!=NULL) { cptr+=5; *cptr='B'; }
    else if ((cptr=strstr(spec->ss.abfile,"e2ds_B"))!=NULL) { cptr+=5; *cptr='A'; }
    else {
      warnmsg("UVES_r2Dspec_harps(): File name for spectrum\n\t%s\n\
\tdoes not contain string 'e2ds_[AB]' so sky spectrum file\n\
\tname cannot be determined. Better to input the original HARPS\n\
\tdata reduction pipeline products with original names.",spec->file);
      spec->skysub=0; status=0;
    }
    if (spec->skysub) sprintf(spec->ss.file,"%s%s",spec->ss.path,spec->ss.abfile);
  }

  /* Read in blaze filename */
  /* Object flux */
  sprintf(key,"HIERARCH %s DRS BLAZE FILE",observatory); FITS_RKEYS("");
  sprintf(blzfile,"%s%s",spec->path,cards);

  /* Determine the arm we're using, the nominal central wavelength,
     the slit-width, CCD binning, spatial pixel scale, the temperature
     and the atmospheric pressure */
  /* At the moment, some of these values are just set to arbitrary numbers */
  spec->cwl=spec->sw=spec->pixscal=-1.0;
  spec->binx=spec->biny=-1;
  /* This may not be the temperature at the grating, which is the T value we really should record */
  if (!spec->fvers) {
    sprintf(key,"HIERARCH %s INS TEMP23 VAL",observatory); FITS_RKEYD("");
    spec->temp=cardd;
    sprintf(key,"HIERARCH %s TEL AMBI PRES START",observatory); FITS_RKEYD("");
    spec->pres=cardd;
    sprintf(key,"HIERARCH %s TEL AMBI PRES END",observatory); FITS_RKEYD("");
    spec->pres+=cardd; spec->pres*=0.5;
  } else {
    sprintf(key,"HIERARCH %s INS VVOUTCOL_T MEAN",observatory); FITS_RKEYD("");
    spec->temp=cardd;
    sprintf(key,"HIERARCH %s METEO PRESSURE",observatory); FITS_RKEYD("");
    spec->pres=cardd;
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
	errormsg("UVES_r2Dspec_harps(): Cannot read format of keyword\n\
%s='%s' in FITS file\n\t%s",key,date,spec->file);
      /* Convert hours, mins and secs into the UT */
      spec->ut=sec/3600.0+min/60.0+hour;
    } else errormsg("UVES_r2Dspec_harps(): Cannot read format of keyword\n\
%s='%s' in FITS file\n\t%s",key,date,spec->file);

    /* Get latitude, longitude and altitude of observatory */
    if (!spec->fvers) {
      FITS_RKEYD("HIERARCH ESO TEL GEOLAT"); spec->lat=cardd;
      FITS_RKEYD("HIERARCH ESO TEL GEOLON"); spec->lon=cardd;
      /* Convert Longitude to normal convention. ESO files use East as
	 positive longitudes. However, the normal convention is for West
	 to be positive. */
      if (spec->lon<0.0) spec->lon*=-1.0;
      else spec->lon=360.0-spec->lon;
      FITS_RKEYD("HIERARCH ESO TEL GEOELEV"); spec->alt=cardd;
    } else {
      FITS_RKEYS("GEOLAT"); 
      if (sscanf(cards,"%lf %lf %lf %s",&deg,&min,&sec,hem)!=4)
	errormsg("UVES_r2Dspec_harps(): Cannot interpret header card\n\
\t%s='%s' from FITS file\n\t%s.",key,cards,spec->file);
      spec->lat=sec/3600.0+min/60.0+deg;
      if (!strncmp(hem,"S",1) || !strncmp(hem,"s",1)) spec->lat*=-1.0;
      FITS_RKEYS("GEOLON"); 
      if (sscanf(cards,"%lf %lf %lf %s",&deg,&min,&sec,hem)!=4)
	errormsg("UVES_r2Dspec_harps(): Cannot interpret header card\n\
\t%s='%s' from FITS file\n\t%s.",key,cards,spec->file);
      spec->lon=sec/3600.0+min/60.0+deg;
      /* Convert Longitude to normal convention, i.e. west is positive. */
      if (!strncmp(hem,"E",1) || !strncmp(hem,"e",1)) spec->lon*=-1.0;
      else { FITS_RKEYD("GEOELEV"); }
      spec->alt=cardd;
    }

    /* Get object RA and DEC and equinox */
    if (!spec->fvers) { FITS_RKEYD("RA"); }
    else { FITS_RKEYD("RA-DEG"); }
    spec->ra=cardd;
    /* Convert RA into hours */
    spec->ra/=15.0;
    if (!spec->fvers) { FITS_RKEYD("DEC"); }
    else { FITS_RKEYD("DEC-DEG"); }
    spec->dec=cardd;
    FITS_RKEYD("EQUINOX"); spec->equ=cardd;
  }

  /* Get exposure time */
  FITS_RKEYD("EXPTIME"); spec->etime=cardd;

  /* Get e-/ADU conversion ratio and read noise per pixel */
  sprintf(key,"HIERARCH %s DRS CCD CONAD",observatory); FITS_RKEYD("");
  eperADU=cardd;
  sprintf(key,"HIERARCH %s DRS CCD SIGDET",observatory); FITS_RKEYD("");
  rdnoise=cardd;

  /* Get number of echelle orders & allocate memory for echelle order array */
  FITS_RKEYI("NAXIS2"); spec->nor=cardi;
  if (!(spec->or=(echorder *)malloc((size_t)(spec->nor*sizeof(echorder)))))
    errormsg("UVES_r2Dspec_harps(): Could not allocate memory for echelle\n\
\torder array of size %d for spectrum in file\n\t%s",spec->nor,spec->file);

  /* Find number of pixels to read in for each order */
  FITS_RKEYI("NAXIS1"); spec->or[0].np=cardi;
  for (i=1; i<spec->nor; i++) spec->or[i].np=spec->or[0].np;

  /* Allocate memory for data arrays and fill wavelength, flux arrays */
  if (!UVES_memspec(spec,par,spec->nor,0))
    errormsg("UVES_r2Dspec_harps(): Error returned from UVES_memspec() when\n\
\tattempting to allocate memory for data arrays for file\n\t%s",spec->file);

  /* Get image dimensions */
  if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
    errormsg("UVES_r2Dspec_harps(): Couldn't get image dimensions for\n\
\tfile %s",spec->file);
  if (!naxis) errormsg("UVES_r2Dspec_harps(): Couldn't find input image extension\n\
\tfor file %s",spec->file);
  npixobj=naxes[0];

  /* Read in flux information */
  for (i=0,first=1; i<spec->nor; i++,first+=npixobj) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].fl,&anynul,&status))
      errormsg("UVES_r2Dspec_harps(): Cannot read flux array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }

  /** Attempt to read the wavelength solution information **/
  /* Get the polynomial degree */
  sprintf(key,"HIERARCH %s DRS CAL TH DEG LL",observatory); FITS_RKEYI("");
  ndeg=cardi+1;
  /* Loop over orders */
  for (i=0; i<spec->nor; i++) {
    spec->or[i].nwpol=ndeg;
    /* For each order, read in the wavelength solution polynomial coefficients */
    for (j=0,k=ndeg*i; j<ndeg; j++,k++) {
      sprintf(key,"HIERARCH %s DRS CAL TH COEFF LL%d",observatory,k); FITS_RKEYD("");
      spec->or[i].wpol[j]=cardd;
    }
  }

  /* Close flux fits file */
  fits_close_file(infits,&status);


  /** Open and read in sky image, if available **/
  if (spec->skysub) {
    /* Open input file as FITS file */
    if (fits_open_file(&infits,UVES_replace_envinstr(spec->ss.file),READONLY,&status)) {
      warnmsg("UVES_r2Dspec_harps(): Cannot open FITS file\n\t%s\n\
\twhich nominally corresponds to sky spectrum for spectrum in file\n\
\t%s.",spec->ss.file,spec->file);
      spec->skysub=0; status=0;
    } else {
      /* Check HDU type */
      if (fits_get_hdu_type(infits,&hdutype,&status))
	errormsg("UVES_r2Dspec_harps(): Cannot get HDU type for file\n\t%s",spec->ss.file);
      if (hdutype!=IMAGE_HDU)
	errormsg("UVES_r2Dspec_harps(): File not a FITS image: %s",spec->ss.file);
      /* Check number of HDUs */
      if (fits_get_num_hdus(infits,&hdunum,&status))
	errormsg("UVES_r2Dspec_harps(): Cannot find number of HDUs in file\n\
\t%s",spec->ss.file);
      if (hdunum>1) errormsg("UVES_r2Dspec_harps(): Number of HDUs is %d instead\n\
\tof %d (at most) in file\n\t%s",hdunum,1,spec->ss.file);
      /* Read in blaze filename */
      sprintf(key,"HIERARCH %s DRS BLAZE FILE",observatory); FITS_RKEYS("");
      sprintf(sblzfile,"%s%s",spec->path,cards);
      /* Get number of echelle orders & allocate memory for echelle order array */
      FITS_RKEYI("NAXIS2"); spec->ss.nor=cardi;
      /* Note: It's often the case that the sky spectrum does not
	 contain the same number of extracted echelle orders as the
	 corresponding object spectrum. We're not assuming they're
	 equal in what follows. This implies that UVES_skysub.c has to
	 match echelle orders in the sky and object spectra using the
	 wavelength ranges as a guide */
      if (!(spec->ss.or=(echorder *)malloc((size_t)(spec->ss.nor*sizeof(echorder)))))
	errormsg("UVES_r2Dspec_harps(): Could not allocate memory for echelle\n\
\torder array of size %d for sky spectrum in file\n\t%s",spec->ss.nor,spec->ss.file);
      /* Find number of pixels to read in for each order */
      FITS_RKEYI("NAXIS1"); spec->ss.or[0].np=cardi;
      for (i=1; i<spec->ss.nor; i++) spec->ss.or[i].np=spec->ss.or[0].np;
      /* Allocate memory for data arrays and fill wavelength, flux arrays */
      if (!UVES_memspec(spec,par,spec->ss.nor,1))
	errormsg("UVES_r2Dspec_harps(): Error returned from UVES_memspec() when\n\
\tattempting to allocate memory for data arrays for file\n\t%s",spec->ss.file);
      /* Get image dimensions */
      if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
	errormsg("UVES_r2Dspec_harps(): Couldn't get image dimensions for\n\
\tfourth HDU in file\n\t%s",spec->ss.file);
      if (!naxis) errormsg("UVES_r2Dspec_harps(): Couldn't find input image extension\n\
\tfor file %s",spec->ss.file);
      npixsky=naxes[0];
      /* Read in sky flux information */
      for (i=0,first=1; i<spec->ss.nor; i++,first+=npixsky) {
	if (fits_read_img(infits,TDOUBLE,first,spec->ss.or[i].np,&nulval,
			  spec->ss.or[i].fl,&anynul,&status))
	  errormsg("UVES_r2Dspec_harps(): Cannot read flux array for order\n\
\t%d in file\n\t%s",i+1,spec->ss.file);
      }
      /** Attempt to read the wavelength solution information **/
      /* Get the polynomial degree */
      sprintf(key,"HIERARCH %s DRS CAL TH DEG LL",observatory); FITS_RKEYI("");
      ndeg=cardi+1;
      /* Loop over orders */
      for (i=0; i<spec->ss.nor; i++) {
	spec->ss.or[i].nwpol=ndeg;
	/* For each order, read in the wavelength solution polynomial coefficients */
	for (j=0,k=ndeg*i; j<ndeg; j++,k++) {
	  sprintf(key,"HIERARCH %s DRS CAL TH COEFF LL%d",observatory,k); FITS_RKEYD("");
	  spec->ss.or[i].wpol[j]=cardd;
	}
      }
      /* Close sky FITS file */
      fits_close_file(infits,&status);
    }
  }

  /* Generate error array from flux array */
  /* It is entirely beyond me why the HARPS/HARPN reduction pipelines
     do not construct error arrays and include these in the flux
     files. Nevertheless, by assuming a fixed number of effective
     spatial pixels which are combined to form a single extracted
     spectral pixel (here, 5 pixels), that provides a total read noise
     and so we can re-construct an approximate error array from the
     flux (photon count) array */
  /* TBD: Use a Poisson distribution for object flux. Also, error
     estimate is underestimated for low flux levels, particularly
     those that fluctuate below zero, because the error is being
     estimated directly from flux; use some local smoothing as minimum
     flux level?? */
  /* Guess at effective number of spatial pixels per extracted spectral pixel */
  neffspatpix=(!spec->fvers) ? HARPSNSPATPIX : HARPNNSPATPIX;
  /* Construct error array for the object flux */
  for (i=0; i<spec->nor; i++) {
    for (j=0; j<spec->or[i].np; j++) {
      spec->or[i].er[j]=
	(spec->or[i].fl[j]*eperADU+neffspatpix*rdnoise*rdnoise)/eperADU;
      spec->or[i].er[j]=(spec->or[i].er[j]>0.0) ? sqrt(spec->or[i].er[j]) : 0.0;
    }
  }
  /* Construct error array for the sky flux, if available */
  if (spec->skysub) {
    for (i=0; i<spec->ss.nor; i++) {
      for (j=0; j<spec->ss.or[i].np; j++) {
	spec->ss.or[i].er[j]=
	  (spec->ss.or[i].fl[j]*eperADU+neffspatpix*rdnoise*rdnoise)/eperADU;
	spec->ss.or[i].er[j]=(spec->ss.or[i].er[j]>0.0) ?
	  sqrt(spec->ss.or[i].er[j]) : 0.0;
      }
    }
  }

  /* If ThAr information is to be read in, do it now */
  if (par->thar==1) {
    errormsg("UVES_r2Dspec_harps(): Currently there is no functionality\n\
\tto read in the ThAr spectrum from HARPS files. Sorry.");
  }

  /** Open and read in blaze image and blaze-correct the flux and
      error arrays **/
  /** Object flux **/
  /* Open input file as FITS file */
  if (fits_open_file(&infits,UVES_replace_envinstr(blzfile),READONLY,&status)) {
    warnmsg("UVES_r2Dspec_harps(): Cannot open FITS file\n\t%s\n\
\twhich is nominated as the corresponding blaze file\n\
\tfor the spectrum in file\n\t%s.",blzfile,spec->file);
    status=0;
  } else {
    /* Check HDU type */
    if (fits_get_hdu_type(infits,&hdutype,&status))
      errormsg("UVES_r2Dspec_harps(): Cannot get HDU type for file\n\t%s",blzfile);
    if (hdutype!=IMAGE_HDU)
      errormsg("UVES_r2Dspec_harps(): File not a FITS image: %s",blzfile);
    /* Check number of HDUs */
    if (fits_get_num_hdus(infits,&hdunum,&status))
      errormsg("UVES_r2Dspec_harps(): Cannot find number of HDUs in file\n\
\t%s",blzfile);
    if (hdunum>1) errormsg("UVES_r2Dspec_harps(): Number of HDUs is %d instead\n\
\tof %d (at most) in file\n\t%s",hdunum,1,blzfile);
    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec_harps(): Couldn't get image dimensions for\n\
\tfourth HDU in file\n\t%s",blzfile);
    /* Check that image dimensions are the same as for the flux image */
    if (naxes[0]!=npixobj || naxes[1]!=spec->nor)
      errormsg("UVES_r2Dspec_harps(): Blaze image in file\n\t%s\n\
\thas different dimensions to corresponding spectrum in file\n\t%s",blzfile,spec->file);
    /* Allocate memory for matrix to hold blaze */
    if ((blz=dmatrix(spec->nor,spec->or[0].np))==NULL)
      errormsg("UVES_r2Dspec_harps(): Cannot allocate memory for blz\n\
\tmatrix of size %dx%d for file\n\t%s",spec->nor,spec->or[0].np,blzfile);
    /* Read in blaze information */
    for (i=0,first=1; i<spec->nor; i++,first+=npixobj) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,blz[i],&anynul,
			&status))
	errormsg("UVES_r2Dspec_harps(): Cannot read blaze array for order\n\
\t%d in file\n\t%s",i+1,blzfile);
    }
    /* Close blaze FITS file */
    fits_close_file(infits,&status);
    /* Blaze-correct the flux and error arrays appropriately */
    for (i=0; i<spec->nor; i++) {
      for (j=0; j<spec->or[i].np; j++) {
	if (blz[i][j]>0.0) {
	  spec->or[i].fl[j]/=blz[i][j]; spec->or[i].er[j]/=blz[i][j];
	} else spec->or[i].er[j]=0.0;
      }
    }
    /* Free blaze matrix */
    free(*blz); free(blz);
  }
  /* Sky flux */
  if (spec->skysub) {
    /* Open input file as FITS file */
    if (fits_open_file(&infits,UVES_replace_envinstr(sblzfile),READONLY,&status)) {
      warnmsg("UVES_r2Dspec_harps(): Cannot open FITS file\n\t%s\n\
\twhich is nominated as the corresponding blaze file\n\
\tfor the sky spectrum in file\n\t%s.\n\
\ttWill not sky-subtract the corresponding object spectrum.",sblzfile,spec->ss.file);
      spec->skysub=0; status=0;
    } else {
      /* Check HDU type */
      if (fits_get_hdu_type(infits,&hdutype,&status))
	errormsg("UVES_r2Dspec_harps(): Cannot get HDU type for file\n\t%s",sblzfile);
      if (hdutype!=IMAGE_HDU)
	errormsg("UVES_r2Dspec_harps(): File not a FITS image: %s",sblzfile);
      /* Check number of HDUs */
      if (fits_get_num_hdus(infits,&hdunum,&status))
	errormsg("UVES_r2Dspec_harps(): Cannot find number of HDUs in file\n\
\t%s",sblzfile);
      if (hdunum>1) errormsg("UVES_r2Dspec_harps(): Number of HDUs is %d instead\n\
\tof %d (at most) in file\n\t%s",hdunum,1,sblzfile);
      /* Get image dimensions */
      if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
	errormsg("UVES_r2Dspec_harps(): Couldn't get image dimensions for\n\
\tfourth HDU in file\n\t%s",sblzfile);
      /* Check that image dimensions are the same as for the sky flux image */
      if (naxes[0]!=npixsky || naxes[1]!=spec->ss.nor)
	errormsg("UVES_r2Dspec_harps(): Blaze image in file\n\t%s\n\
\thas different dimensions to corresponding spectrum in file\n\t%s",
		 sblzfile,spec->ss.file);
      /* Allocate memory for matrix to hold blaze */
      if ((blz=dmatrix(spec->ss.nor,spec->ss.or[0].np))==NULL)
	errormsg("UVES_r2Dspec_harps(): Cannot allocate memory for blz\n\
\tmatrix of size %dx%d for file\n\t%s",spec->ss.nor,spec->ss.or[0].np,sblzfile);
      /* Read in blaze information */
      for (i=0,first=1; i<spec->ss.nor; i++,first+=npixsky) {
	if (fits_read_img(infits,TDOUBLE,first,spec->ss.or[i].np,&nulval,blz[i],&anynul,
			  &status))
	  errormsg("UVES_r2Dspec_harps(): Cannot read blaze array for order\n\
\t%d in file\n\t%s",i+1,sblzfile);
      }
      /* Close blaze FITS file */
      fits_close_file(infits,&status);
      /* Blaze-correct the sky flux and error arrays appropriately */
      for (i=0; i<spec->ss.nor; i++) {
	for (j=0; j<spec->ss.or[i].np; j++) {
	  if (blz[i][j]>0.0) {
	    spec->ss.or[i].fl[j]/=blz[i][j]; spec->ss.or[i].er[j]/=blz[i][j];
	  } else spec->ss.or[i].er[j]=0.0;
	}
      }
      /* Free blaze matrix */
      free(*blz); free(blz);
    }
  }


  /* Set order identity numbers and parameters in flux spectrum */
  for (i=0; i<spec->nor; i++) {
    spec->or[i].id=i; spec->or[i].sid=spec->id;
    /* Find the number of useful pixels in the order */
    if (par->thar<=1) {
      for (spec->or[i].nuse=0,j=0; j<spec->or[i].np; j++) {
	if (spec->or[i].er[j]>0.0) spec->or[i].nuse++;
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
	errormsg("UVES_r2Dspec_harps(): Error returned from UVES_memspec() when\n\
\tattempting to free memory from order %d of file\n\t%s",i+1,spec->file);
    }
  }

  /* Set order identity numbers and parameters in sky spectrum */
  if (spec->skysub) {
    for (i=0; i<spec->ss.nor; i++) {
      spec->ss.or[i].id=i; spec->ss.or[i].sid=spec->id;
      /* Find the number of useful pixels in the order */
      for (spec->ss.or[i].nuse=0,j=0; j<spec->ss.or[i].np; j++) {
	if (spec->ss.or[i].er[j]>0.0) spec->ss.or[i].nuse++;
	else spec->ss.or[i].st[j]=RCLIP;
      }
      if (spec->ss.or[i].nuse<MINUSE) {
	/* Free the allocated memory for the raw arrays */
	if (!UVES_memspec(spec,par,-i-1,1))
	  errormsg("UVES_r2Dspec_harps(): Error returned from UVES_memspec() when\n\
\tattempting to free memory from order %d of file\n\t%s",i+1,spec->ss.file);
      }
    }
  }

  /* Make sure this spectrum will be combined into combined spectrum later */
  spec->comb=1;

  return 1;

}
