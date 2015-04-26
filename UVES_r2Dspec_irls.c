/****************************************************************************
* Read in a FITS long-slit spectrum that has been reduced with IRAF
****************************************************************************/

#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <longnam.h>
#include "UVES_popler.h"
#include "astron.h"
#include "memory.h"
#include "error.h"

int UVES_r2Dspec_irls(spectrum *spec, params *par) {

  double   hour=0.0,min=0.0,sec=0.0;
  double   nulval=0.0;
  long     naxes[9] = {0,0,0,0,0,0,0,0,0};
  int      npix=0;
  int      hdutype=0,hdunum=0,status=0,bitpix=0,first=1,naxis=0,anynul=0;
  int      i=0,j=0,k=0;
  int      idum=0;
  char     comment[FLEN_COMMENT]="\0";
  char     date[FLEN_KEYWORD]="\0";
  char     time[FLEN_KEYWORD]="\0";
  char     *cptr;
  fitsfile *infits;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,spec->file,READONLY,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  if (fits_get_hdu_type(infits,&hdutype,&status))
    errormsg("UVES_r2Dspec_irls(): Cannot get HDU type for file\n\t%s",spec->file);
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_irls(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec_irls(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (hdunum!=1)
    errormsg("UVES_r2Dspec_irls(): Number of HDUs is %d instead of %d\n\
\tin file %s",hdunum,1,spec->file);

  /* Get object name */
  if (par->thar<=1) {
    if (fits_read_key(infits,TSTRING,"OBJECT",spec->obj,comment,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read value of header\n\
\tcard %s from FITS file\n\t%s.","OBJECT",spec->file);
  }
  else sprintf(spec->obj,"thar_wav");
  /* Alter object name to remove some special characters */
  while ((cptr=strchr(spec->obj,' '))!=NULL) *cptr='_';
  while ((cptr=strchr(spec->obj,'+'))!=NULL) *cptr='p';
  while ((cptr=strchr(spec->obj,'-'))!=NULL) *cptr='m';

  /* Read in entire main header as array of strings */
  if (fits_get_hdrspace(infits,&(spec->nhead_ori),NULL,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read number of keys in main\n\
\theader in FITS file\n\t%s.",spec->file);
  if ((spec->head_ori=cmatrix(spec->nhead_ori,FLEN_CARD))==NULL)
    errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for head_ori\n\
\tmatrix of size %dx%d for file %s\n",spec->nhead_ori,FLEN_CARD,spec->file);
  for (i=0; i<spec->nhead_ori; i++) {
    if (fits_read_record(infits,i+1,spec->head_ori[i],&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read header key %d from\n\
\tmain header of file %s\n",i+1,spec->file);
  }

  /* Read in or, in this case, set archival filename */
  sprintf(spec->arfile,"%s",spec->abfile);

  /* Determine the arm we're using, the nominal central wavelength,
     the slit-width, CCD binning, spatial pixel scale, the temperature
     and the atmospheric pressure */
  /* At the moment, these values are just set to arbitrary numbers */
  spec->cwl=spec->sw=spec->pixscal=spec->temp=spec->pres=-1.0;
  spec->binx=spec->biny=-1;

  if (par->thar<=1) {
    /* Get modified julian day */
    if (fits_read_key(infits,TDOUBLE,"MJD-OBS",&(spec->jd),comment,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->file);
    /* Convert to Julian day */
    spec->jd+=2400000.5;
    /* Get date of observation and convert to year, month and day */
    if (fits_read_key(infits,TSTRING,"DATE-OBS",date,comment,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DATE-OBS",spec->file);
    if (strlen(date)>10) {
      /* This means that the UT is appended to the date */
      if (sscanf(date,"%d-%d-%dT%lf:%lf:%lf",&(spec->year),&(spec->month),
		 &(spec->day),&hour,&min,&sec)!=6)
	errormsg("UVES_r2Dspec_irls(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->file);
      /* Convert hours, mins and secs into the UT */
      spec->ut=sec/3600.0+min/60.0+hour;
    } else {
      /* Only the date should be given */
      if (sscanf(date,"%d-%d-%d",&(spec->year),&(spec->month),&(spec->day))!=3)
	errormsg("UVES_r2Dspec_irls(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->file);
      /* Get universal time */ 
      if (fits_read_key(infits,TSTRING,"UTSTART",time,comment,&status))
	errormsg("UVES_r2Dspec_irls(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","UTSTART",spec->file);
      /* Convert to hours */
      if (sscanf(time,"%lf:%lf:%lf",&hour,&min,&sec)!=3)
	errormsg("UVES_r2Dspec_irls(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","UTSTART",time,spec->file);
      spec->ut=sec/3600.0+min/60.0+hour;
    }

    /* Get latitude, longitude and altitude of observatory */
    if (fits_read_key(infits,TDOUBLE,"LAT_OBS",&(spec->lat),comment,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","LAT_OBS",spec->file);
    /* Convert Longitude to normal convention. It is assumed that the
       longitude in the header is in the usual astronomical (and IRAF)
       convention of being positive for logitudes West of Greenwich */
    if (fits_read_key(infits,TDOUBLE,"LONG_OBS",&(spec->lon),comment,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","LONG_OBS",spec->file);
    if (fits_read_key(infits,TDOUBLE,"ALT_OBS",&(spec->alt),comment,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","ALT_OBS",spec->file);
    /** Get object RA and DEC and equinox **/
    if (fits_read_key(infits,TDOUBLE,"RA",&(spec->ra),comment,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","RA",spec->file);
    /* Convert RA into hours */
    spec->ra/=15.0;
    if (fits_read_key(infits,TDOUBLE,"DEC",&(spec->dec),comment,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DEC",spec->file);
    if (fits_read_key(infits,TDOUBLE,"EQUINOX",&(spec->equ),comment,&status)) {
      warnmsg("UVES_r2Dspec_irls(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.\n\tAssuming equinox = %6.1lf","EQUINOX",spec->file,
	      C_J2000); spec->equ=C_J2000; status=0;
    }
  }

  /* Get exposure time */
  if (fits_read_key(infits,TDOUBLE,"EXPTIME",&(spec->etime),comment,&status))
    errormsg("UVES_r2Dspec_irls(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","EXPTIME",spec->file);

  /* Assume that there is one 'order' per file and allocate memory for
     echelle order array */
  spec->nor=1;
  if (!(spec->or=(echorder *)malloc((size_t)(spec->nor*sizeof(echorder)))))
    errormsg("Could not allocate memory for echelle order array of size %d",
	     spec->nor);

  /* Find number of pixels to read in for each order */
  if (fits_read_key(infits,TINT,"NAXIS1",&(spec->or[0].np),comment,&status))
    errormsg("UVES_r2Dspec_irls(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS1",spec->file);
  for (i=1; i<spec->nor; i++) spec->or[i].np=spec->or[0].np;

  /* Allocate memory for data arrays and fill wavelength, flux arrays */
  for (i=0; i<spec->nor; i++) {
    if ((spec->or[i].vhwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for vac. wavel.\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].vhrwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for right-edge\n\
\tvac-wavel. array of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,
	       spec->file);
    if ((spec->or[i].fl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for flux\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].er=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for error\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].res=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for resolution\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if (par->thar==1) {
      if ((spec->or[i].th=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for ThAr\n\
\tflux array of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
      if ((spec->or[i].ter=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for ThAr\n\
\terror array of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,
		 spec->file);
    }
    if ((spec->or[i].st=iarray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for status\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if (par->thar==1) {
      if ((spec->or[i].tst=iarray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_irls(): Cannot allocate memory for ThAr status\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    }

    /* Initialise status array */
    for (j=0; j<spec->or[i].np; j++) spec->or[i].st[j]=1;
    if (par->thar==1) for (j=0; j<spec->or[i].np; j++) spec->or[i].tst[j]=1;
  }

  /* Get image dimensions */
  if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
    errormsg("UVES_r2Dspec_irls(): Couldn't get image dimensions for\n\
\tfile %s",spec->file);
  if (!naxis)
    errormsg("UVES_r2Dspec_irls(): Couldn't find input image extension\n\
\tfor file %s",spec->file);
  npix=naxes[0];

  /* Read in flux information */
  for (i=0,first=1; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].fl,&anynul,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read flux array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }

  /* Close flux fits file */
  fits_close_file(infits,&status);

  /* Attempt to open error flux fits file */
  if (par->thar<=1) {
    if (fits_open_file(&infits,spec->erfile,READONLY,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot open FITS file\n\t%s",spec->erfile);

    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec_irls(): Couldn't get image dimensions for\n\
\tfile %s",spec->erfile);
    if (!naxis)
      errormsg("UVES_r2Dspec_irls(): Couldn't find input image extension\n\
\tfor file %s",spec->erfile);

    /* Some weak checks that this really is the right error array */
    if (naxes[0]!=npix)
      errormsg("UVES_r2Dspec_irls(): Mismatch in array sizes between\n\
\tflux and error files\n\t%s &\n\t%s",spec->file,spec->erfile);
    if (naxes[1]!=0 && naxes[1]!=1)
      errormsg("UVES_r2Dspec_irls(): Mismatch in number of orders between\n\
\tflux and error files\n\t%s &\n\t%s",spec->file,spec->erfile);

    /* Read in error information */
    for (i=0,first=1; i<spec->nor; i++,first+=npix) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
			spec->or[i].er,&anynul,&status))
	errormsg("UVES_r2Dspec_irls(): Cannot read error array for order\n\
\t%d in file\n\t%s",i+1,spec->erfile);
    }

    /* Close error flux fits file */
    fits_close_file(infits,&status);
  }

  /* Scale flux and error up by 10^18 to avoid rounding error intolerance */
  for (i=0; i<spec->nor; i++) {
    for (j=0; j<spec->or[i].np; j++) {
      spec->or[i].fl[j]*=1.e18;
      spec->or[i].er[j]*=1.e18;
    }
  }

  /* If ThAr information is to be read in, do it now */
  if (par->thar==1) {
    if (fits_open_file(&infits,spec->thfile,READONLY,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot open FITS file\n\t%s",spec->thfile);

    /* Read in or, in this case, set archival filename */
    sprintf(spec->tharfile,"%s",spec->abthfile);

    /* Get modified julian day */
    if (fits_read_key(infits,TDOUBLE,"MJD-OBS",&(spec->thjd),comment,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->thfile);
    /* Convert to Julian day */
    spec->thjd+=2400000.5;

    /* Find the temperature in each arm and the atmospheric pressure */
    /* At the moment, these values are just set to arbitrary numbers */
    spec->thtemp=spec->thpres=-1.0;

    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec_irls(): Couldn't get image dimensions for\n\
\tfile %s",spec->thfile);
    if (!naxis) errormsg("UVES_r2Dspec_irls(): Couldn't find input image\n\
\textension for file\n\t%s",spec->thfile);

    /* Some weak checks that this really is the right ThAr array */
    if (naxes[0]!=npix)
      errormsg("UVES_r2Dspec_irls(): Mismatch in array sizes between\n\
\tflux and error files\n\t%s &\n\t%s",spec->file,spec->thfile);
    if (naxes[1]!=spec->nor)
      errormsg("UVES_r2Dspec_irls(): Mismatch in number of orders between\n\
\tflux and error files\n\t%s &\n\t%s",spec->file,spec->thfile);

    /* Read in ThAr information */
    for (i=0,first=1; i<spec->nor; i++,first+=npix) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
			spec->or[i].th,&anynul,&status))
	errormsg("UVES_r2Dspec_irls(): Cannot read ThAr flux array for order\n\
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
	if (spec->or[i].er[j]>DRNDTOL) spec->or[i].nuse++;
	else spec->or[i].st[j]=RCLIP;
	if (par->thar==1) {
	  if (spec->or[i].th[j]==0.0) spec->or[i].tst[j]=RCLIP;
	}
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
      free(spec->or[i].vhwl); free(spec->or[i].vhrwl); free(spec->or[i].fl);
      free(spec->or[i].er); free(spec->or[i].res); free(spec->or[i].st);
    }
      if (par->thar==1) {
	free(spec->or[i].th); free(spec->or[i].ter); free(spec->or[i].tst);
      }
  }

  /* Attempt to open wavelength solution fits file */
  if (fits_open_file(&infits,spec->wlfile,READONLY,&status))
      errormsg("UVES_r2Dspec_irls(): Cannot open FITS file\n\t%s",spec->wlfile);
  
  /* Read in wavelength solution information */
  if (fits_read_key(infits,TINT,"DC-FLAG",&idum,comment,&status))
    errormsg("UVES_r2Dspec_irls(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DC-FLAG",spec->file);
  if (idum) errormsg("UVES_r2Dspec_irls(): File %s does not contain a linear\n\
\twavelength solution, required when using IRAF LSS data.",spec->file);

  /* Read in linear wavelength solution to wpol parameters in terms of
     polynomial soluton */
  for (i=0; i<spec->nor; i++) {
    if (fits_read_key(infits,TDOUBLE,"CRVAL1",&(spec->or[i].wpol[0]),comment,
		      &status))
    errormsg("UVES_r2Dspec_irls(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","CRVAL1",spec->file);
    if (fits_read_key(infits,TDOUBLE,"CD1_1",&(spec->or[i].wpol[1]),comment,
		      &status))
    errormsg("UVES_r2Dspec_irls(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","CD1_1",spec->file);
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
