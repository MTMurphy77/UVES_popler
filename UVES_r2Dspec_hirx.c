/****************************************************************************
* Read in a 2D FITS echelle spectrum that has been reduced with HIRES_REDUX
****************************************************************************/

#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <longnam.h>
#include "UVES_popler.h"
#include "fit.h"
#include "memory.h"
#include "error.h"
#include "const.h"

int UVES_r2Dspec_hirx(spectrum *spec, params *par) {

  double   hour=0.0,min=0.0,sec=0.0,onorm=0.0;
  double   nulval=0.0;
  double   nrmt[2];
  double   dwl=0.0;
  double   *blzfit=NULL;
  double   *coeff=NULL,*coeffo=NULL;
  double   **blz=NULL;
  long     nrows=0,No=0;
  long     naxes[9]={0,0,0,0,0,0,0,0,0};
  int      npix=0,col=0,norm=0;
  int      hdutype=0,hdunum=0,status=0,bitpix=0,first=1,naxis=0,anynul=0;
  int      i=0,j=0,k=0;
  int      *ord=NULL;
  char     comment[FLEN_COMMENT]="\0";
  char     date[FLEN_KEYWORD]="\0",time[FLEN_KEYWORD]="\0";
  char     *cptr;
  fitsfile *infits;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,UVES_replace_envinstr(spec->file),READONLY,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  if (fits_get_hdu_type(infits,&hdutype,&status))
    errormsg("UVES_r2Dspec_hirx(): Cannot get HDU type for file\n\t%s",spec->file);
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_hirx(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec_hirx(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (par->thar<=1) {
    if (hdunum==4) norm=1;
    else if (hdunum==5) norm=0;
    else errormsg("UVES_r2Dspec_hirx(): Number of HDUs is %d but it\n\
\tmust be either %d or %d in file\n\t%s",hdunum,4,5,spec->file);
  }

  /* Get object name */
  if (par->thar<=1) {
    if (fits_read_key(infits,TSTRING,"OBJECT",spec->obj,comment,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","OBJECT",spec->file);
  } else sprintf(spec->obj,"thar_wav");
  /* Alter object name to remove some special characters */
  while ((cptr=strchr(spec->obj,' '))!=NULL) *cptr='_';
  while ((cptr=strchr(spec->obj,'+'))!=NULL) *cptr='p';
  while ((cptr=strchr(spec->obj,'-'))!=NULL) *cptr='m';

  /* Read in entire main header as array of strings */
  if (fits_get_hdrspace(infits,&(spec->nhead_ori),NULL,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read number of keys in main\n\
\theader in FITS file\n\t%s.",spec->file);
  if ((spec->head_ori=cmatrix(spec->nhead_ori,FLEN_CARD))==NULL)
    errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for head_ori\n\
\tmatrix of size %dx%d for file %s\n",spec->nhead_ori,FLEN_CARD,spec->file);
  for (i=0; i<spec->nhead_ori; i++) {
    if (fits_read_record(infits,i+1,spec->head_ori[i],&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read header key %d from\n\
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
    if (fits_read_key(infits,TDOUBLE,"MJD",&(spec->jd),comment,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD",spec->file);
    /* Convert to Julian day */
    spec->jd+=2400000.5;
    /* Get date of observation and convert to year, month and day */
    if (fits_read_key(infits,TSTRING,"DATE-OBS",date,comment,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DATE-OBS",spec->file);
    if (strlen(date)>10) {
      /* This means that the UT is appended to the date */
      if (sscanf(date,"%d-%d-%dT%lf:%lf:%lf",&(spec->year),&(spec->month),
		 &(spec->day),&hour,&min,&sec)!=6)
	errormsg("UVES_r2Dspec_hirx(): Cannot read format of keyword\n\
%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->file);
      /* Convert hours, mins and secs into the UT */
      spec->ut=sec/3600.0+min/60.0+hour;
    } else {
      /* Only the date should be given */
      if (sscanf(date,"%d-%d-%d",&(spec->year),&(spec->month),&(spec->day))!=3)
	errormsg("UVES_r2Dspec_hirx(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->file);
      /* Get universal time from UTC */ 
      if (fits_read_key(infits,TSTRING,"UTC",time,comment,&status))
	errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","UTC",spec->file);
      /* Convert to hours */
      if (!UVES_hex2dec(time,&(spec->ut)))
	errormsg("UVES_r2Dspec_hirx(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","UTC",time,spec->file);
    }
    /* Get latitude, longitude and altitude of observatory */
    /* At the moment, these are hard-coded assuming that we really are
       dealing with Keck/HIRES data */
    spec->lat=19.8259; spec->lon=155.4751; spec->alt=4160.0;
    /* Get object RA and DEC and equinox */
    if (fits_read_key(infits,TSTRING,"RA",time,comment,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","RA",spec->file);
    /* Convert RA into hours */
    if (!UVES_hex2dec(time,&(spec->ra)))
      errormsg("UVES_r2Dspec_hirx(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","RA",time,spec->file);
    if (fits_read_key(infits,TSTRING,"DEC",time,comment,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DEC",spec->file);
    /* Convert to degrees */
    if (!UVES_hex2dec(time,&(spec->dec)))
      errormsg("UVES_r2Dspec_hirx(): Cannot read format of keyword\n\
%s=%s in FITS file\n\t%s","DEC",time,spec->file);
    if (fits_read_key(infits,TDOUBLE,"EQUINOX",&(spec->equ),comment,&status))
      warnmsg("UVES_r2Dspec_hirx(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.\n\tAssuming equinox = 2000.0","EQUINOX",spec->file);
  }

  /* Get exposure time */
  if (fits_read_key(infits,TDOUBLE,"EXPTIME",&(spec->etime),comment,&status))
    errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","EXPTIME",spec->file);

  /* Read heliocentric velicity correction */
  /* Note: HIRES REDUX files are assumed to have already been
     corrected to the vacuum-heliocentric frame. However, if
     atmospheric lines are to be masked out by the user, then this
     needs to be done in the observer's frame, so the masked regions
     need to have the same heliocentric correction applied to them as
     the input spectra have */
  if (fits_read_key(infits,TDOUBLE,"HELVEL",&(spec->vhel),comment,&status))
    errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","HELVEL",spec->file);

  /* Get number of echelle orders & allocate memory for echelle order array */
  if (fits_read_key(infits,TINT,"NAXIS2",&(spec->nor),comment,&status))
    errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS2",spec->file);
  if (!(spec->or=(echorder *)malloc((size_t)(spec->nor*sizeof(echorder)))))
    errormsg("Could not allocate memory for echelle order array of size %d",
	     spec->nor);

  /* Find number of pixels to read in for each order */
  if (fits_read_key(infits,TINT,"NAXIS1",&(spec->or[0].np),comment,&status))
    errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS1",spec->file);
  for (i=1; i<spec->nor; i++) spec->or[i].np=spec->or[0].np;

  /* Allocate memory for data arrays and fill wavelength, flux arrays */
  for (i=0; i<spec->nor; i++) {
    /* Unlike most other instruments/pipelines, MAGE files have
       an actual wavelength array to be read in, so allocate memory
       for the raw wavelength array first */
    if ((spec->or[i].wl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for raw wavel.\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].vhwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for vac. wavel.\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].vhrwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for right-edge\n\
\tvac-wavel. array of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,
	       spec->file);
    if ((spec->or[i].fl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for flux\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].er=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for error\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].res=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for resolution\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if (par->thar==1) {
      if ((spec->or[i].th=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for ThAr\n\
\tflux array of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
      if ((spec->or[i].ter=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for ThAr\n\
\terror array of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,
		 spec->file);
    }
    if ((spec->or[i].st=iarray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for status\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if (par->thar==1) {
      if ((spec->or[i].tst=iarray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for ThAr status\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    }

    /* Initialise status array */
    for (j=0; j<spec->or[i].np; j++) spec->or[i].st[j]=1;
    if (par->thar==1) for (j=0; j<spec->or[i].np; j++) spec->or[i].tst[j]=1;
  }

  /* Get image dimensions */
  if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
    errormsg("UVES_r2Dspec_hirx(): Couldn't get image dimensions for\n\
\tfile %s",spec->file);
  if (!naxis) errormsg("UVES_r2Dspec_hirx(): Couldn't find input image extension\n\
\tfor file %s",spec->file);
  npix=naxes[0];

  /* Read in flux information */
  for (i=0,first=1; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].fl,&anynul,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read flux array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }

  /** Next HDU assumed to contain error information **/
  if (par->thar<=1) {
    /* Move to next HDU */
    if (fits_movrel_hdu(infits,1,&hdutype,&status))
      errormsg("UVES_r2Dspec_hirx(): Could not move to second HDU\n\
\tin file %s",spec->file);
    /* Check HDU type */
    if (hdutype!=IMAGE_HDU)
      errormsg("UVES_r2Dspec_hirx(): Second extension not a FITS image\n\
\tin file\n\t%s",spec->file);
    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec_hirx(): Couldn't get image dimensions for\n\
\tsecond HDU in file\n\t%s",spec->file);
    if (!naxis) errormsg("UVES_r2Dspec_hirx(): Couldn't find input image extension\n\
\tfor second HDU in file\n\t%s",spec->file);
    /* Check that image dimensions are the same as for the flux image */
    if (naxes[0]!=npix || naxes[1]!=spec->nor)
      errormsg("UVES_r2Dspec_hirx(): Second extension has different\n\
\tdimensions to first in file\n\t%s",spec->file);
    /* Read in error information */
    for (i=0,first=1; i<spec->nor; i++,first+=npix) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,spec->or[i].er,
			&anynul,&status))
	errormsg("UVES_r2Dspec_hirx(): Cannot read error array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
      /* The HIRES_REDUX format uses variance, not sigma, so take sqrt here */
      for (j=0; j<spec->or[i].np; j++)
	spec->or[i].er[j]=(spec->or[i].er[j]>0.0) ? sqrt(spec->or[i].er[j]) : -INFIN;
    }
  }

  /** Next HDU assumed to contain wavelength information **/
  if (par->thar<=1) {
    /* Move to next HDU */
    if (fits_movrel_hdu(infits,1,&hdutype,&status))
      errormsg("UVES_r2Dspec_hirx(): Could not move to third HDU\n\
\tin file %s",spec->file);
    /* Check HDU type */
    if (hdutype!=IMAGE_HDU)
      errormsg("UVES_r2Dspec_hirx(): Third extension not a FITS image\n\
\tin file\n\t%s",spec->file);
    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec_hirx(): Couldn't get image dimensions for\n\
\tthird HDU in file\n\t%s",spec->file);
    if (!naxis) errormsg("UVES_r2Dspec_hirx(): Couldn't find input image extension\n\
\tfor third HDU in file\n\t%s",spec->file);
    /* Check that image dimensions are the same as for the flux image */
    if (naxes[0]!=npix || naxes[1]!=spec->nor)
      errormsg("UVES_r2Dspec_hirx(): Third extension has different\n\
\tdimensions to first in file\n\t%s",spec->file);
    /* Read in raw wavelength information */
    for (i=0,first=1; i<spec->nor; i++,first+=npix) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,spec->or[i].wl,
			&anynul,&status))
	errormsg("UVES_r2Dspec_hirx(): Cannot read wavelength array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
    }
  }
  /* Must treat the special case for HIRX files where wavelengths are
     set to zero if there's no information from the chip at those
     wavelengths. At the moment, the fix is just to extrapolate the
     wavelength scale from valid sections of each order out to the
     order edges, only replacing zero wavelength entries */
  for (i=0; i<spec->nor; i++) {
    /* Find first non-zero wavelength element */
    for (j=0; j<spec->or[i].np-1; j++) if (spec->or[i].wl[j]>DRNDTOL) break;
    dwl=spec->or[i].wl[j+1]-spec->or[i].wl[j];
    for (k=j-1; k>=0; k--) spec->or[i].wl[k]=spec->or[i].wl[k+1]-dwl;
    /* Find last non-zero wavelength element */
    for (j=spec->or[i].np-1; j>=1; j--) if (spec->or[i].wl[j]>DRNDTOL) break;
    dwl=spec->or[i].wl[j]-spec->or[i].wl[j-1];
    for (k=j+1; k<spec->or[i].np; k++) spec->or[i].wl[k]=spec->or[i].wl[k-1]+dwl;
  }


  /** When there's 5 HDU's, this next HDU should contain the extracted
      flat-field flux (i.e. some measure of the blaze function) **/
  if (par->thar<=1 && !norm) {
    /* Move to next HDU */
    if (fits_movrel_hdu(infits,1,&hdutype,&status))
      errormsg("UVES_r2Dspec_hirx(): Could not move to fourth HDU\n\
\tin file %s",spec->file);
    /* Check HDU type */
    if (hdutype!=IMAGE_HDU)
      errormsg("UVES_r2Dspec_hirx(): Fourth extension not a FITS image\n\
\tin file\n\t%s",spec->file);
    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec_hirx(): Couldn't get image dimensions for\n\
\tfourth HDU in file\n\t%s",spec->file);
    if (!naxis) errormsg("UVES_r2Dspec_hirx(): Couldn't find input image extension\n\
\tfor fourth HDU in file\n\t%s",spec->file);
    /* Check that image dimensions are the same as for the flux image */
    if (naxes[0]!=npix || naxes[1]!=spec->nor)
      errormsg("UVES_r2Dspec_hirx(): Fourth extension has different\n\
\tdimensions to first in file\n\t%s",spec->file);
    /* Allocate memory for matrix to hold blaze */
    if ((blz=dmatrix(spec->nor,spec->or[0].np))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for blz\n\
\tmatrix of size %dx%d for file\n\t%s",spec->nor,spec->or[0].np,spec->file);
    /* Read in blaze information */
    for (i=0,first=1; i<spec->nor; i++,first+=npix) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,blz[i],&anynul,
			&status))
	errormsg("UVES_r2Dspec_hirx(): Cannot read blaze array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
    }
  }

  /* Close flux fits file */
  fits_close_file(infits,&status);

  /* If ThAr information is to be read in, do it now */
  if (par->thar==1) {
    if (fits_open_file(&infits,UVES_replace_envinstr(spec->thfile),READONLY,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot open FITS file\n\t%s",spec->thfile);
    /* Read in or, in this case, set archival filename */
    sprintf(spec->tharfile,"%s",spec->abthfile);
    /* Get modified julian day */
    if (fits_read_key(infits,TDOUBLE,"MJD",&(spec->wc_jd),comment,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD",spec->thfile);
    /* Convert to Julian day */
    spec->wc_jd+=2400000.5;
    /* Find the temperature in each arm and the atmospheric pressure */
    /* At the moment, these values are just set to arbitrary numbers */
    spec->wc_temp=spec->wc_pres=-1.0;
    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec_hirx(): Couldn't get image dimensions for\n\
\tfile %s",spec->thfile);
    if (!naxis) errormsg("UVES_r2Dspec_hirx(): Couldn't find input image\n\
\textension for file\n\t%s",spec->thfile);
    /* Some weak checks that this really is the right ThAr array */
    if (naxes[0]!=npix)
      errormsg("UVES_r2Dspec_hirx(): Mismatch in array sizes between\n\
\tflux and error files\n\t%s &\n\t%s",spec->file,spec->thfile);
    if (naxes[1]!=spec->nor)
      errormsg("UVES_r2Dspec_hirx(): Mismatch in number of orders between\n\
\tflux and error files\n\t%s &\n\t%s",spec->file,spec->thfile);
    /* Read in ThAr information */
    for (i=0,first=1; i<spec->nor; i++,first+=npix) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,spec->or[i].th,
			&anynul,&status))
	errormsg("UVES_r2Dspec_hirx(): Cannot read ThAr flux array for order\n\
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
	/* Have applied an upper threshold for errors here because
	   HIREDUX doesn't apply a negative error bar to regions which
	   are extracted but which are known to contain no flux */
	/* if (spec->or[i].er[j]>DRNDTOL) spec->or[i].nuse++; */
	if (spec->or[i].er[j]>DRNDTOL && spec->or[i].er[j]<500.0) spec->or[i].nuse++;
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

  /** Normalize by fits to the blaze function if necessary **/
  if (par->thar<=1 && !norm) {
    /* Allocate memory for matrix to hold blaze */
    if ((blzfit=darray(spec->or[0].np))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for blzfit\n\
\tarray of size %d for file\n\t%s",spec->or[0].np,spec->file);
    /* Loop over orders and do the fits to the blaze function */
    for (i=0; i<spec->nor; i++) {
      /* Only operate on useful orders */
      if (spec->or[i].nuse>=MINUSE) {
	/* Fit blaze function */
	/*
	if (!UVES_hirx_blzfit(&(spec->or[i]),&(blz[i][0]),
			      10.0,3.0,25,5,0,0,blzfit)) {
	  nferrormsg("UVES_r2Dspec_hirx(): Error returned from UVES_hirx_blzfit()\n\
\twhen attempting to fit blaze function of order %d of spectrum %d,\n\t%s",i+1,
		     spec->id+1,spec->file); return 0;
	}
	*/
	/* Find the running median of the blaze function */
	if (!UVES_hirx_blzfit(&(spec->or[i]),&(blz[i][0]),0.0,3.0,51,2,0,1,blzfit)) {
	  nferrormsg("UVES_r2Dspec_hirx(): Error returned from UVES_hirx_blzfit()\n\
\twhen attempting a running median on blaze function of order %d of\n\
\tspectrum %d,\n\t%s",i+1,spec->id+1,spec->file); return 0;
	}
	/* Divide order by fit to the blaze */
	for (j=0; j<spec->or[i].np; j++) {
	  /*
	  if (i==9) fprintf(stdout,"%lf  %lf  %lf  %lf  %lf  %d\n",spec->or[i].vhwl[j],spec->or[i].fl[j],spec->or[i].er[j],blz[i][j],blzfit[j],spec->or[i].st[j]);
	  */
	  if (spec->or[i].st[j]==1) {
	    spec->or[i].fl[j]/=blzfit[j]; spec->or[i].er[j]/=blzfit[j];
	  }
	}
      }
    }
    /* Clean up */
    free(blzfit); free(*blz); free(blz);
  }

  /* Read in the polynomial wavelength solutions from the extracted ThAr file.
     The ThAr information in this file is extracted differently to the flux and
     so the wavelength polynomials do not apply to the flux and error arrays
     above. Indeed, above we read in the wavelength scale directly, not via
     polynomials like below */
  if (par->thar==2) {
    /* Attempt to open wavelength solution fits file */
    if (fits_open_file(&infits,UVES_replace_envinstr(spec->wlfile),READONLY,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot open FITS file\n\t%s",spec->wlfile);
    /* Check number of HDUs */
    if (fits_get_num_hdus(infits,&hdunum,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot find number of HDUs in file\n\
\t%s",spec->wlfile);
    if (hdunum!=3)
      errormsg("UVES_r2Dspec_hirx(): Number of HDUs is %d instead of %d\n\
\tin file %s",hdunum,spec->nor+1,spec->wlfile);
    /* Move to next HDU */
    if (fits_movrel_hdu(infits,1,&hdutype,&status))
      errormsg("UVES_r2Dspec_hirx(): Could not move to second HDU\n\
\tin file %s",spec->wlfile);
    /* Check HDU type */
    if (hdutype!=BINARY_TBL)
      errormsg("UVES_r2Dspec_hirx(): Extension %d not a binary table\n\
\tin file\n\t%s",2,spec->wlfile);
    /* Chebyshev coefficients for each echelle order need to be
       constructed from 2D fit coefficients */
    /* Find number of axes to be read in */
    if (fits_get_num_cols(infits,&naxis,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read number of columns\n\
\t%s in FITS file\n\t%s.","TFIELDS",spec->wlfile);
    if (naxis!=5) errormsg("UVES_r2Dspec_hirx(): The binary table in file\n\t%s\n\
\thas %d columns. It should have 5",spec->wlfile,naxis);
    /* Find number of rows to be read in */
    if (fits_get_num_rows(infits,&nrows,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read number of rows\n\
\t%s in FITS file\n\t%s.","NAXIS2",spec->wlfile);
    if (nrows!=1) errormsg("UVES_r2Dspec_hirx(): The binary table in file\n\t%s\n\
\thas %d rows. It should have 1",spec->wlfile,nrows);
    /* Find and read in the normalizing coefficients
       for the dispersion direction */
    if (fits_get_colnum(infits,CASEINSEN,"NRM",&col,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","NRM",spec->wlfile);
    if (fits_read_col(infits,TDOUBLE,col,1,1,2,&nulval,&(spec->or[0].wpol[NWPOL-2]),
		      &anynul,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read column '%s'\n\
\tin binary table in FITS file\n\t%s","NRM",spec->wlfile);
    /* Write these values into last elements of all echelle orders' wpol arrays */
    for (i=1; i<spec->nor; i++) {
      spec->or[i].wpol[NWPOL-2]=spec->or[0].wpol[NWPOL-2];
      spec->or[i].wpol[NWPOL-1]=spec->or[0].wpol[NWPOL-1];
    }
    /* Find and read in the normalizing coefficients for the spatial
       (order) dierction */
    if (fits_get_colnum(infits,CASEINSEN,"NRMT",&col,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","NRMT",spec->wlfile);
    if (fits_read_col(infits,TDOUBLE,col,1,1,2,&nulval,nrmt,&anynul,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read column '%s'\n\
\tin binary table in FITS file\n\t%s","NRMT",spec->wlfile);
    /* Find and read in the order of the Legendre polynomial fit in the
       dispersion direction */
    if (fits_get_colnum(infits,CASEINSEN,"NY",&col,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","NY",spec->wlfile);
    if (fits_read_col(infits,TINT,col,1,1,1,&nulval,&(spec->or[0].nwpol),&anynul,
		      &status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read column '%s'\n\
\tin binary table in FITS file\n\t%s","NY",spec->wlfile);
    /* Write this value into all echelle orders */
    for (i=1; i<spec->nor; i++) spec->or[i].nwpol=spec->or[0].nwpol;
    /* Find and read in the order of the Legendre polynomial fit in the
       spatial direction */
    if (fits_get_colnum(infits,CASEINSEN,"NO",&col,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","NO",spec->wlfile);
    if (fits_read_col(infits,TINT,col,1,1,1,&nulval,&No,&anynul,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read column '%s'\n\
\tin binary table in FITS file\n\t%s","NO",spec->wlfile);
    /* Make sure wpol arrays are big enough to handle this number of
       coefficients */
    if (NWPOL<spec->or[0].nwpol+2)
      errormsg("UVES_r2Dspec_hirx(): Order of wavelength calibration\n\
\tpolynomial in spectral direction (=%d) must be <= NWPOL-2=%d for\n	\
\tHIRES REDUX data. Try increasing NWPOL to %d and recompiling",spec->or[0].nwpol,
	       NWPOL-2,spec->or[0].nwpol+2);
    /* Find and read in the coefficient matrix */
    if (fits_get_colnum(infits,CASEINSEN,"RES",&col,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","RES",spec->wlfile);
    if (fits_get_coltype(infits,col,&i,&(naxes[0]),&(naxes[1]),&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot find type of column\n\
\tnamed '%s' in binary table in FITS file\n\t%s","RES",spec->wlfile);
    if (naxes[0]!=spec->or[0].nwpol*No)
      errormsg("UVES_r2Dspec_hirx(): Size of coefficient matrix (=%d)\n\
\tdoes not equal polynomial order in dispersion x spatial\n	       \
\tdirection, %dx%d, in binary table in FITS file\n\t%s",naxes[0],
	       spec->or[0].nwpol,No,spec->wlfile);
    /* Allocate memory for coefficient array */
    if ((coeff=darray(naxes[0]))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for coeff\n\
\tarray of size %d",naxes[0]+1);
    if (fits_read_col(infits,TDOUBLE,col,1,1,naxes[0],&nulval,coeff,&anynul,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read column '%s'\n\
\tin binary table in FITS file\n\t%s","RES",spec->wlfile);
    /* Diffraction order numbers contained in the next HDU */
    /* Move to next HDU */
    if (fits_movrel_hdu(infits,1,&hdutype,&status))
      errormsg("UVES_r2Dspec_hirx(): Could not move to third HDU\n\
\tin file %s",spec->wlfile);
    /* Check HDU type */
    if (hdutype!=BINARY_TBL)
      errormsg("UVES_r2Dspec_hirx(): Extension %d not a binary table\n\
\tin file\n\t%s",i+1,spec->wlfile);
    /* Find number of axes to be read in */
    if (fits_get_num_cols(infits,&naxis,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read number of columns\n\
\t%s in FITS file\n\t%s.","NAXIS2",spec->wlfile);
    /* Find number of rows to be read in */
    if (fits_get_num_rows(infits,&nrows,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read number of rows\n\
\t%s in FITS file\n\t%s.","NAXIS2",spec->wlfile);
    if (nrows!=spec->nor)
      errormsg("UVES_r2Dspec_hirx(): Second binary table in file\n\t%s\n\
\thas %d rows. It should be equal to the number of echelle orders, %d\n	\
\tin this case",spec->wlfile,nrows,spec->nor);
    /* Allocate memory for order number array */
    if ((ord=iarray(nrows))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for ord\n\
\tarray of size %d",nrows);
    /* Find and read in the diffraction order numbers */
    if (fits_get_colnum(infits,CASEINSEN,"ORDER",&col,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","ORDER",spec->wlfile);
    if (fits_read_col(infits,TINT,col,1,1,nrows,&nulval,ord,&anynul,&status))
      errormsg("UVES_r2Dspec_hirx(): Cannot read column '%s'\n\
\tin binary table in FITS file\n\t%s","ORDER",spec->wlfile);
    /* Allocate memory for Legendre coefficients */
    if ((coeffo=darray(No))==NULL)
      errormsg("UVES_r2Dspec_hirx(): Cannot allocate memory for coeffo\n\
\tarray of size %d",No);
    /* Generate Legendre coefficients for single orders from 2D fit coefficients */
    for (i=0; i<spec->nor; i++) {
      /* Normalize order number */
      onorm=2.0*((double)ord[i]-nrmt[0])/nrmt[1];
      if (!svdfit_legendre(onorm,coeffo,No))
	errormsg("UVES_r2Dspec_hirx(): Error returned from svdfit_legendre()\n\
\twhen calculating coefficients for order %d of file\n\t%s",i+1,spec->wlfile);
      /* Loop to generate the j'th Legendre coefficient for this order */
      for (j=0; j<spec->or[i].nwpol; j++) {
	/* Sum of products of coefficient matrix elements and Legendre series */
	for (k=No-1,spec->or[i].wpol[j]=0.0; k>=0; k--)
	  spec->or[i].wpol[j]+=coeff[j*No+k]*coeffo[k];
	spec->or[i].wpol[j]/=(double)ord[i];
      }
    }
    /* Clean up */
    free(coeff); free(coeffo); free(ord);
    /* Close wavelength solution fits file */
    fits_close_file(infits,&status);
  }

  /* Set median resolution to zero for lack of any other information */
  spec->arcfwhm=0.0;

  /* Make sure this spectrum will be combined into combined spectrum later */
  spec->comb=1;

  return 1;

}
