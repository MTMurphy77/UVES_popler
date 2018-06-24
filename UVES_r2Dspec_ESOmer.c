/****************************************************************************
* Read in a 2D FITS echelle spectrum that has been reduced with IRAF
****************************************************************************/

#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <longnam.h>
#include "UVES_popler.h"
#include "astron.h"
#include "memory.h"
#include "error.h"

int UVES_r2Dspec_ESOmer(spectrum *spec, params *par) {

  double   hour=0.0,min=0.0,sec=0.0;
  double   nulval=0.0;
  long     naxes[9] = {0,0,0,0,0,0,0,0,0};
  int      npix=0;
  int      hdutype=0,hdunum=0,status=0,bitpix=0,first=1,naxis=0,anynul=0;
  int      i=0,j=0,k=0;
  char     card[FLEN_CARD]="\0";
  char     comment[FLEN_COMMENT]="\0";
  char     date[FLEN_KEYWORD]="\0";
  char     *cptr;
  fitsfile *infits;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,spec->file,READONLY,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  if (fits_get_hdu_type(infits,&hdutype,&status))
    errormsg("UVES_r2Dspec_ESOmer(): Cannot get HDU type for file\n\t%s",spec->file);
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_ESOmer(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec_ESOmer(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (hdunum!=1)
    errormsg("UVES_r2Dspec_ESOmer(): Number of HDUs is %d instead of %d\n\
\tin file %s",hdunum,1,spec->file);

  /* Read in entire main header as array of strings */
  if (fits_get_hdrspace(infits,&(spec->nhead_ori),NULL,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read number of keys in main\n\
\theader in FITS file\n\t%s.",spec->file);
  if ((spec->head_ori=cmatrix(spec->nhead_ori,FLEN_CARD))==NULL)
    errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for head_ori\n\
\tmatrix of size %dx%d for file %s\n",spec->nhead_ori,FLEN_CARD,spec->file);
  for (i=0; i<spec->nhead_ori; i++) {
    if (fits_read_record(infits,i+1,spec->head_ori[i],&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read header key %d from\n\
\tmain header of file %s\n",i+1,spec->file);
  }

  /* Get object name */
  if (par->thar<=1) {
    if (fits_read_key(infits,TSTRING,"HIERARCH ESO OBS TARG NAME",spec->obj,comment,
		      &status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header\n\
\tcard %s from FITS file\n\t%s.","HIERARCH ESO OBS TARG NAME",spec->file);
  }
  else sprintf(spec->obj,"thar_wav");
  /* Alter object name to remove some special characters */
  while ((cptr=strchr(spec->obj,' '))!=NULL) *cptr='_';
  while ((cptr=strchr(spec->obj,'+'))!=NULL) *cptr='p';
  while ((cptr=strchr(spec->obj,'-'))!=NULL) *cptr='m';

  /* Read in archival filename */
  if (fits_read_key(infits,TSTRING,"ARCFILE",spec->arfile,comment,&status))
    errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","ARCFILE",spec->file);

  /* Determine the arm we're using, the nominal central wavelength,
     the slit-width, CCD binning, spatial pixel scale, the temperature
     and the atmospheric pressure */
  if (fits_read_key(infits,TSTRING,"HIERARCH ESO INS PATH",card,comment,&status))
    errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO INS PATH",spec->file);
  if (strstr(card,"RED")!=NULL || strstr(card,"red")!=NULL) {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS TEMP2 MEAN",&(spec->temp),
		      comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS TEMP2 MEAN",spec->file);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS GRAT2 WLEN",&(spec->cwl),
		      comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS GRAT2 WLEN",spec->file);
  } else {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS TEMP1 MEAN",&(spec->temp),
		      comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS TEMP1 MEAN",spec->file);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS GRAT1 WLEN",&(spec->cwl),
		      comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS GRAT1 WLEN",spec->file);
  }
  if (spec->cwl<UVESBORR) {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SLIT2 WID",&(spec->sw),comment,
		      &status))
	errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SLIT2 WID",spec->file);
  } else {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SLIT3 WID",&(spec->sw),comment,
		      &status))
	errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SLIT3 WID",spec->file);
  }
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SENS26 MEAN",&(spec->pres),
		    comment,&status))
    errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SENS26 MEAN",spec->file);

  if (par->thar<=1) {
    /* Get modified julian day */
    if (fits_read_key(infits,TDOUBLE,"MJD-OBS",&(spec->jd),comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->file);
    /* Convert to Julian day */
    spec->jd+=2400000.5;
    /* Get date of observation and convert to year, month and day */
    if (fits_read_key(infits,TSTRING,"DATE-OBS",date,comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DATE-OBS",spec->file);
    if (strlen(date)>10) {
      /* This means that the UT is appended to the date */
      if (sscanf(date,"%d-%d-%dT%lf:%lf:%lf",&(spec->year),&(spec->month),
		 &(spec->day),&hour,&min,&sec)!=6)
	errormsg("UVES_r2Dspec_ESOmer(): Cannot read format of keyword\n\
%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->file);
      /* Convert hours, mins and secs into the UT */
      spec->ut=sec/3600.0+min/60.0+hour;
    } else {
      /* Only the date should be given */
      if (sscanf(date,"%d-%d-%d",&(spec->year),&(spec->month),&(spec->day))!=3)
	errormsg("UVES_r2Dspec_ESOmer(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->file);
      /* Get universal time from TM-START or UTC */ 
      if (fits_read_key(infits,TDOUBLE,"TM-START",&(spec->ut),comment,&status)) {
	if (fits_read_key(infits,TDOUBLE,"UTC",&(spec->ut),comment,&status))
	  errormsg("UVES_r2Dspec_ESOmer(): Cannot read values of header cards\n\
\t%s & %s from FITS file\n\t%s.","TM-START","UTC",spec->file);
      }
      /* Convert to hours */
      spec->ut/=3600.0;
    }
    /* Get latitude, longitude and altitude of observatory */
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL GEOLAT",&(spec->lat),
		      comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO TEL GEOLAT",spec->file);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL GEOLON",&(spec->lon),
		      comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO TEL GEOLON",spec->file);
    /* Convert Longitude to normal convention. UVES files use East as
       positive longitudes. However, the normal convention is for West
       to be positive. */
    if (spec->lon<0.0) spec->lon*=-1.0;
    else spec->lon=360.0-spec->lon;
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL GEOELEV",&(spec->alt),
		      comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO TEL GEOELEV",spec->file);
    /* Get object RA and DEC and equinox */
    if (fits_read_key(infits,TDOUBLE,"RA",&(spec->ra),comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","RA",spec->file);
    /* Convert RA into hours */
    spec->ra/=15.0;
    if (fits_read_key(infits,TDOUBLE,"DEC",&(spec->dec),comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DEC",spec->file);
    if (fits_read_key(infits,TDOUBLE,"EQUINOX",&(spec->equ),comment,&status)) {
      warnmsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.\n\tAssuming equinox = %6.1lf","EQUINOX",spec->file,
	      C_J2000); spec->equ=C_J2000; status=0;
    }
  }

  /* Get exposure time */
  if (fits_read_key(infits,TDOUBLE,"EXPTIME",&(spec->etime),comment,&status))
    errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","EXPTIME",spec->file);

  /* Assume that there is one 'order' per file and allocate memory for
     echelle order array */
  spec->nor=1;
  if (!(spec->or=(echorder *)malloc((size_t)(spec->nor*sizeof(echorder)))))
    errormsg("Could not allocate memory for echelle order array of size %d",
	     spec->nor);

  /* Find number of pixels to read in for each order */
  if (fits_read_key(infits,TINT,"NAXIS1",&(spec->or[0].np),comment,&status))
    errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS1",spec->file);
  for (i=1; i<spec->nor; i++) spec->or[i].np=spec->or[0].np;

  /* Allocate memory for data arrays and fill wavelength, flux arrays */
  for (i=0; i<spec->nor; i++) {
    if ((spec->or[i].vhwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for vac. wavel.\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].vhrwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for right-edge\n\
\tvac-wavel. array of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,
	       spec->file);
    if ((spec->or[i].fl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for flux\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].er=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for error\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].res=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for resolution\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if (par->thar==1) {
      if ((spec->or[i].th=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for ThAr flux\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
      if ((spec->or[i].ter=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for ThAr error\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    }
    if ((spec->or[i].st=iarray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for status\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if (par->thar==1) {
      if ((spec->or[i].tst=iarray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_ESOmer(): Cannot allocate memory for ThAr status\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    }
    /* Initialise status array */
    for (j=0; j<spec->or[i].np; j++) spec->or[i].st[j]=1;
    if (par->thar==1) for (j=0; j<spec->or[i].np; j++) spec->or[i].tst[j]=1;
  }

  /* Get image dimensions */
  if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
    errormsg("UVES_r2Dspec_ESOmer(): Couldn't get image dimensions for\n\
\tfile %s",spec->file);
  if (!naxis) errormsg("UVES_r2Dspec_ESOmer(): Couldn't find input image extension\n\
\tfor file %s",spec->file);
  npix=naxes[0];

  /* Read in flux information */
  for (i=0,first=1; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].fl,&anynul,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read flux array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }

  /* Close flux fits file */
  fits_close_file(infits,&status);

  /* Attempt to open error flux fits file */
  if (par->thar<=1) {
    if (fits_open_file(&infits,spec->erfile,READONLY,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot open FITS file\n\t%s",spec->erfile);

    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Couldn't get image dimensions for\n\
\tfile %s",spec->erfile);
    if (!naxis)
      errormsg("UVES_r2Dspec_ESOmer(): Couldn't find input image extension\n\
\tfor file %s",spec->erfile);

    /* Some weak checks that this really is the right error array */
    if (naxes[0]!=npix)
      errormsg("UVES_r2Dspec_ESOmer(): Mismatch in array sizes between\n\
\tflux and error files\n\t%s &\n\t%s",spec->file,spec->erfile);
    if (naxes[1]!=spec->nor)
      errormsg("UVES_r2Dspec_ESOmer(): Mismatch in number of orders between\n\
\tflux and error files\n\t%s &\n\t%s",spec->file,spec->erfile);

    /* Read in error information */
    for (i=0,first=1; i<spec->nor; i++,first+=npix) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
			spec->or[i].er,&anynul,&status))
	errormsg("UVES_r2Dspec_ESOmer(): Cannot read error array for order\n\
\t%d in file\n\t%s",i+1,spec->erfile);
    }

    /* Close error flux fits file */
    fits_close_file(infits,&status);
  }

  /* If ThAr information is to be read in, do it now */
  if (par->thar==1) {
    if (fits_open_file(&infits,spec->thfile,READONLY,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot open FITS file\n\t%s",spec->thfile);
    /* Read in archival filename */
    if (fits_read_key(infits,TSTRING,"ARCFILE",spec->tharfile,comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","ARCFILE",spec->thfile);
    /* Get modified julian day */
    if (fits_read_key(infits,TDOUBLE,"MJD-OBS",&(spec->wc_jd),comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->thfile);
    /* Convert to Julian day */
    spec->wc_jd+=2400000.5;
    /* Find the temperature in each arm and the atmospheric pressure */
    if (fits_read_key(infits,TSTRING,"HIERARCH ESO INS PATH",card,comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO INS PATH",spec->thfile);
    if (strstr(card,"RED")!=NULL || strstr(card,"red")!=NULL)
      sprintf(card,"HIERARCH ESO INS TEMP2 MEAN");
    else sprintf(card,"HIERARCH ESO INS TEMP1 MEAN");
    if (fits_read_key(infits,TDOUBLE,card,&(spec->wc_temp),comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file %s.",card,spec->thfile);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SENS26 MEAN",
		      &(spec->wc_pres),comment,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SENS26 MEAN",spec->thfile);
    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Couldn't get image dimensions for\n\
\tfile %s",spec->thfile);
    if (!naxis) errormsg("UVES_r2Dspec_ESOmer(): Couldn't find input image extension\n\
\tfor file %s",spec->thfile);
    /* Some weak checks that this really is the right ThAr array */
    if (naxes[0]!=npix)
      errormsg("UVES_r2Dspec_ESOmer(): Mismatch in array sizes between\n\
\tflux and ThAr files\n\t%s &\n\t%s",spec->file,spec->thfile);
    if (naxes[1]!=spec->nor)
      errormsg("UVES_r2Dspec_ESOmer(): Mismatch in number of orders between\n\
\tflux and ThAr files\n\t%s &\n\t%s",spec->file,spec->thfile);
    /* Read in ThAr information */
    for (i=0,first=1; i<spec->nor; i++,first+=npix) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
			spec->or[i].th,&anynul,&status))
	errormsg("UVES_r2Dspec_ESOmer(): Cannot read ThAr flux array for order\n\
\t%d in file\n\t%s",i+1,spec->thfile);
    }
    /* Close ThAr flux fits file */
    fits_close_file(infits,&status);
  }

  /* Set order identity numbers and parameters */
  for (i=0; i<spec->nor; i++) {
    spec->or[i].id=i; spec->or[i].sid=spec->id;
    /* Find the number of useful pixels in the order */
    if (par->thar<=1) {
      for (spec->or[i].nuse=0,j=0; j<spec->or[i].np; j++) {
	if (spec->or[i].er[j]>DRNDTOL && spec->or[i].er[j]!=1.0) spec->or[i].nuse++;
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
      free(spec->or[i].vhwl); free(spec->or[i].vhrwl); free(spec->or[i].fl);
      free(spec->or[i].er); free(spec->or[i].res); free(spec->or[i].st);
      if (par->thar==1) {
	free(spec->or[i].th); free(spec->or[i].ter); free(spec->or[i].tst);
      }
    }
  }

  /* Attempt to open wavelength solution fits file */
  if (fits_open_file(&infits,spec->wlfile,READONLY,&status))
      errormsg("UVES_r2Dspec_ESOmer(): Cannot open FITS file\n\t%s",spec->wlfile);
  
  /* Read in linear wavelength solution to wpol parameters in terms of
     polynomial soluton */
  for (i=0; i<spec->nor; i++) {
    if (fits_read_key(infits,TDOUBLE,"CRVAL1",&(spec->or[i].wpol[0]),comment,
		      &status))
    errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","CRVAL1",spec->file);
    if (fits_read_key(infits,TDOUBLE,"CDELT1",&(spec->or[i].wpol[1]),comment,
		      &status))
    errormsg("UVES_r2Dspec_ESOmer(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","CDELT1",spec->file);
    spec->or[i].nwpol=2;
  }
  /* Close wavelength solution fits file */
  fits_close_file(infits,&status);

  /* Set median resolution to zero for lack of any other information */
  spec->arcfwhm=0.0;

  /* Make sure this spectrum will be combined into combined spectrum later */
  spec->comb=1;

  return 1;

}
