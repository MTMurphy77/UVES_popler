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

int UVES_r2Dspec(spectrum *spec, params *par) {

  double   hour=0.0,min=0.0,sec=0.0,mnorm=0.0,dm=0.0;
  double   nulval=0.0;
  double   *coefd=NULL,*resol=NULL,*resid=NULL,*seeing=NULL;
  long     nrows=0;
  long     naxes[9]={0,0,0,0,0,0,0,0,0};
  int      npix=0,col=0,ncoefd=0,No=0,ms=0,me=0,mslope=0,idum=0;
  int      hdutype=0,hdunum=0,hdumove=0,status=0,bitpix=0,first=1,naxis=0,anynul=0;
  int      i=0,j=0,k=0,l=0,m=0;
  int      *ord1=NULL,*ord2=NULL;
  char     comment[FLEN_COMMENT]="\0";
  char     date[FLEN_KEYWORD]="\0";
  char     card[FLEN_CARD]="\0",dummy[FLEN_KEYWORD]="\0";
  char     *inclist[1],search[FLEN_CARD]="HISTORY\0";
  char     *cptr;
  statset  stat;
  fitsfile *infits;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,spec->file,READONLY,&status))
      errormsg("UVES_r2Dspec(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  if (fits_get_hdu_type(infits,&hdutype,&status))
    errormsg("UVES_r2Dspec(): Cannot get HDU type for file\n\t%s",spec->file);
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (hdunum>2)
    errormsg("UVES_r2Dspec(): Number of HDUs is %d instead\n\
\tof %d (at most) in file\n\t%s",hdunum,2,spec->file);

  /* Get file version to indicate whether this is a MIDAS- or CPL-produced file */
  fits_read_key(infits,TSTRING,"ORIGIN",card,comment,&status);
  if (!status && (strstr(card,"MIDAS")!=NULL || strstr(card,"midas")!=NULL))
    spec->fvers=0;
  else {
    /* For now, first advance one HDU when attempting to deal with CPL-written files */
    status=0;    
    if (fits_movrel_hdu(infits,1,&hdutype,&status))
      errormsg("UVES_r2Dspec(): Cannot move to second HDU\n\
\tin FITS file\n\t%s.",spec->file);
    if (fits_read_key(infits,TSTRING,"HIERARCH ESO PRO REC1 DRS ID",card,comment,
		      &status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","HIERARCH ESO PRO REC1 DRS ID",spec->file);
    if (strstr(card,"CPL")!=NULL || strstr(card,"cpl")!=NULL) spec->fvers=1;
    else errormsg("UVES_r2Dspec(): No indication that file\n\t%s\n\
\toriginates from MIDAS reduction from header keyword %s\n\
\tbut no indication of CPL reduction from header keyword\n\t%s",
		  spec->file,"ORIGIN","HIERARCH ESO PRO REC1 DRS ID");
  }

  /* Read in entire main header as array of strings */
  if (fits_get_hdrspace(infits,&(spec->nhead_ori),NULL,&status))
      errormsg("UVES_r2Dspec(): Cannot read number of keys in main\n\
\theader in FITS file\n\t%s.",spec->file);
  if ((spec->head_ori=cmatrix(spec->nhead_ori,FLEN_CARD))==NULL)
    errormsg("UVES_r2Dspec(): Cannot allocate memory for head_ori\n\
\tmatrix of size %dx%d for file %s\n",spec->nhead_ori,FLEN_CARD,spec->file);
  for (i=0; i<spec->nhead_ori; i++) {
    if (fits_read_record(infits,i+1,spec->head_ori[i],&status))
      errormsg("UVES_r2Dspec(): Cannot read header key %d from\n\
\tmain header of file %s\n",i+1,spec->file);
  }

  /* Get object name */
  if (par->thar<=1) {
    if (fits_read_key(infits,TSTRING,"HIERARCH ESO OBS TARG NAME",spec->obj,
		      comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","HIERARCH ESO OBS TARG NAME",spec->file);
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
  if (!spec->fvers &&
      fits_read_key(infits,TSTRING,"ARCFILE",spec->arfile,comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","ARCFILE",spec->file);
  else if (spec->fvers==1) {
    if (fits_read_key(infits,TSTRING,"MJD-OBS",card,comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->file);
    if ((cptr=strchr(comment,'('))==NULL)
      errormsg("UVES_r2Dspec(): Cannot find start of\n\
\tUT time stamp in comment of header card\n\t%s in file\n\t%s","MJD-OBS",spec->file);
    sprintf(spec->arfile,"UVES.%s",cptr+1);
    if ((cptr=strchr(spec->arfile,')'))==NULL)
      errormsg("UVES_r2Dspec(): Cannot find end of\n\
\tUT time stamp in comment of header card\n\t%s in file\n\t%s","MJD-OBS",spec->file);
    sprintf(cptr,"%s",".fits");
  }

  /* Determine the arm we're using, the instrument & derotator mode,
     the nominal central wavelength, the slit-width, CCD binning,
     spatial pixel scale, airmass, seeing, the temperature and the
     atmospheric pressure */
  if (fits_read_key(infits,TSTRING,"HIERARCH ESO INS PATH",spec->inspath,comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO INS PATH",spec->file);
  if (fits_read_key(infits,TSTRING,"HIERARCH ESO INS MODE",spec->insmode,comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO INS MODE",spec->file);
  if (fits_read_key(infits,TSTRING,"HIERARCH ESO INS DROT MODE",spec->insdrot,comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO INS DROT MODE",spec->file);
  if (strstr(spec->inspath,"RED")!=NULL || strstr(spec->inspath,"red")!=NULL) {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS TEMP2 MEAN",&(spec->temp),
		      comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS TEMP2 MEAN",spec->file);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS GRAT2 WLEN",&(spec->cwl),
		      comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS GRAT2 WLEN",spec->file);
  } else {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS TEMP1 MEAN",&(spec->temp),
		      comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS TEMP1 MEAN",spec->file);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS GRAT1 WLEN",&(spec->cwl),
		      comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS GRAT1 WLEN",spec->file);
  }
  if (spec->cwl<UVESBORR) {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SLIT2 WID",&(spec->sw),comment,
		      &status))
	errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SLIT2 WID",spec->file);
    if (fits_read_key(infits,TINT,"HIERARCH ESO INS GRAT1 ENC",&(spec->encoder),
		      comment,&status)) {
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS GRAT1 ENC",spec->file);
    }
  } else {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SLIT3 WID",&(spec->sw),comment,
		      &status))
	errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SLIT3 WID",spec->file);
    if (fits_read_key(infits,TINT,"HIERARCH ESO INS GRAT2 ENC",&(spec->encoder),
		      comment,&status)) {
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS GRAT2 ENC",spec->file);
    }
  }
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS PIXSCALE",&(spec->pixscal),
		    comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO INS PIXSCALE",spec->file);
  if (fits_read_key(infits,TINT,"HIERARCH ESO DET WIN1 BINY",&(spec->binx),comment,
		    &status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO DET WIN1 BINY",spec->file);
  if (fits_read_key(infits,TINT,"HIERARCH ESO DET WIN1 BINX",&(spec->biny),comment,
		    &status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO DET WIN1 BINX",spec->file);
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL AMBI FWHM START",&(spec->seeing[0]),
		    comment,&status)) {
    spec->seeing[0]=-1.0; status=0.0;
  }
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL AMBI FWHM END",&(spec->seeing[1]),
		    comment,&status)) {
    spec->seeing[1]=-1.0; status=0.0;
  }
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL AIRM START",&(spec->airmass[0]),
		    comment,&status)) {
    spec->airmass[1]=-1.0; status=0;
  }
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL AIRM END",&(spec->airmass[1]),
		    comment,&status)) {
    spec->airmass[1]=-1.0; status=0;
  }
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SENS26 MEAN",&(spec->pres),
		    comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SENS26 MEAN",spec->file);
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO QC VRAD HELICOR",&(spec->vhel_head),
		    comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO QC VRAD HELICOR",spec->file);
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO QC VRAD BARYCOR",&(spec->vbar_head),
		    comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO QC VRAD BARYCOR",spec->file);
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO GEN MOON RA",&(spec->moon_ra),
		    comment,&status)) {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL MOON RA",&(spec->moon_ra),
		      comment,&status)) { spec->moon_ra=-1.0; status=0.0; }
  }
  if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO GEN MOON DEC",&(spec->moon_dec),
		    comment,&status)) {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL MOON DEC",&(spec->moon_dec),
		      comment,&status)) { spec->moon_dec=-91.0; status=0.0; }
  } else spec->moon_dec/=15.0;
  if (spec->moon_ra>0.0 && spec->moon_dec>-90.0) {
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO GEN MOON DIST",
		      &(spec->moon_ang),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO GEN MOON DIST",spec->file);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO GEN MOON PHASE",
		      &(spec->moon_phase),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO GEN MOON PHASE",spec->file);
  } else { spec->moon_ang=-1.0; spec->moon_phase=-1.0; }

  if (par->thar<=1) {
    /* Get modified julian day */
    if (fits_read_key(infits,TDOUBLE,"MJD-OBS",&(spec->jd),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->file);
    /* Convert to Julian day */
    spec->jd+=2400000.5;
    /* Get date of observation and convert to year, month and day */
    if (fits_read_key(infits,TSTRING,"DATE-OBS",date,comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DATE-OBS",spec->file);
    if (strlen(date)>10) {
      /* This means that the UT is appended to the date */
      if (sscanf(date,"%d-%d-%dT%lf:%lf:%lf",&(spec->year),&(spec->month),
		 &(spec->day),&hour,&min,&sec)!=6)
	errormsg("UVES_r2Dspec(): Cannot read format of keyword\n\
%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->file);
      /* Convert hours, mins and secs into the UT */
      spec->ut=sec/3600.0+min/60.0+hour;
    } else {
      /* Only the date should be given */
      if (sscanf(date,"%d-%d-%d",&(spec->year),&(spec->month),&(spec->day))!=3)
	errormsg("UVES_r2Dspec(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->file);
      /* Get universal time from TM-START or UTC */ 
      if (fits_read_key(infits,TDOUBLE,"TM-START",&(spec->ut),comment,&status)) {
	if (fits_read_key(infits,TDOUBLE,"UTC",&(spec->ut),comment,&status))
	  errormsg("UVES_r2Dspec(): Cannot read values of header cards\n\
\t%s & %s from FITS file\n\t%s.","TM-START","UTC",spec->file);
      }
      /* Convert to hours */
      spec->ut/=3600.0;
    }
    /* Get latitude, longitude and altitude of observatory */
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL GEOLAT",&(spec->lat),
		      comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO TEL GEOLAT",spec->file);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL GEOLON",&(spec->lon),
		      comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO TEL GEOLON",spec->file);
    /* Convert Longitude to normal convention. UVES files use East as
       positive longitudes. However, the normal convention is for West
       to be positive. */
    if (spec->lon<0.0) spec->lon*=-1.0;
    else spec->lon=360.0-spec->lon;
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO TEL GEOELEV",&(spec->alt),
		      comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO TEL GEOELEV",spec->file);
    /* Get object RA and DEC and equinox */
    if (fits_read_key(infits,TDOUBLE,"RA",&(spec->ra),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","RA",spec->file);
    /* Convert RA into hours */
    spec->ra/=15.0;
    if (fits_read_key(infits,TDOUBLE,"DEC",&(spec->dec),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DEC",spec->file);
    if (fits_read_key(infits,TDOUBLE,"EQUINOX",&(spec->equ),comment,&status)) {
      warnmsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.\n\tAssuming equinox = %6.1lf","EQUINOX",spec->file,
	      C_J2000); spec->equ=C_J2000; status=0;
    }
  }

  /* Get exposure time and encoder step */
  if (fits_read_key(infits,TDOUBLE,"EXPTIME",&(spec->etime),comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","EXPTIME",spec->file);

  /* For now we have to move back to first HDU if dealing with CPL-written files */
  if (spec->fvers==1 && fits_movabs_hdu(infits,1,&hdutype,&status))
    errormsg("UVES_r2Dspec(): Cannot move to first HDU\n\
\tin FITS file\n\t%s.",spec->file);
  /* Read in some basic header card information to define program & OB IDs, type of exposure.
   heliocentric and barycentric corrections */
  if (fits_read_key(infits,TSTRING,"HIERARCH ESO DPR TECH",spec->dprtech,comment,&status)) {
    status=0; sprintf(spec->dprtech,"Not_available");
  }
  if (fits_read_key(infits,TSTRING,"HIERARCH ESO DPR TYPE",spec->dprtype,comment,&status)) {
    status=0; sprintf(spec->dprtype,"Not_available");
  }
  if (fits_read_key(infits,TSTRING,"HIERARCH ESO DPR CATG",spec->dprcatg,comment,&status)) {
    status=0; sprintf(spec->dprcatg,"Not_available");
  }
  /* The ESO OBS ID and PROG ID cards normally sits in the primary HDU but are
     sometimes only present in the first extension in CPL-reduced files */
  if (fits_read_key(infits,TSTRING,"HIERARCH ESO OBS PROG ID",spec->progid,comment,
		    &status)) {
    if (spec->fvers==1) {
      hdumove=1; status=0.0; if (fits_movabs_hdu(infits,2,&hdutype,&status))
	errormsg("UVES_r2Dspec(): Cannot move to second HDU\n\
\tin FITS file\n\t%s.",spec->file);
      if (fits_read_key(infits,TSTRING,"HIERARCH ESO OBS PROG ID",spec->progid,comment,
			&status))
	errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO OBS PROG ID",spec->file);
    }
  }
  if (fits_read_key(infits,TINT,"HIERARCH ESO OBS ID",&(spec->obid),comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO OBS ID",spec->file);
  if (hdumove==1) { if (fits_movabs_hdu(infits,1,&hdutype,&status)) {
		    errormsg("UVES_r2Dspec(): Cannot move to first HDU\n\
\tin FITS file\n\t%s.",spec->file); hdumove=0;
    }
  }

  /* Get number of echelle orders & allocate memory for echelle order array */
  if (fits_read_key(infits,TINT,"NAXIS2",&(spec->nor),comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS2",spec->file);
  if (!(spec->or=(echorder *)malloc((size_t)(spec->nor*sizeof(echorder)))))
    errormsg("UVES_r2Dspec(): Could not allocate memory for echelle\n\
\torder array of size %d",spec->nor);

  /* Find number of pixels to read in for each order */
  if (fits_read_key(infits,TINT,"NAXIS1",&(spec->or[0].np),comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS1",spec->file);
  for (i=1; i<spec->nor; i++) spec->or[i].np=spec->or[0].np;

  /* Allocate memory for data arrays and fill wavelength, flux arrays */
  if (!UVES_memspec(spec,par,spec->nor,0))
    errormsg("UVES_r2Dspec(): Error returned from UVES_memspec() when\n\
\tattempting to allocate memory for data arrays for file\n\t%s",spec->file);

  /* Get image dimensions */
  if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
    errormsg("UVES_r2Dspec(): Couldn't get image dimensions for\n\
\tfile %s",spec->file);
  if (!naxis) errormsg("UVES_r2Dspec(): Couldn't find input image extension\n\
\tfor file %s",spec->file);
  npix=naxes[0];

  /* Read in flux information */
  for (i=0,first=1; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].fl,&anynul,&status))
      errormsg("UVES_r2Dspec(): Cannot read flux array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }

  /** For CPL UVES files, read in the seeing determined from the
      optimal extraction **/
  if (spec->fvers==1 && par->thar<=1) {
    /* Move back to 2nd HDU */
    if (fits_movrel_hdu(infits,1,&hdutype,&status))
      errormsg("UVES_r2Dspec(): Cannot move to second HDU\n\tin FITS file\n\t%s.",
	       spec->file);
    /* Loop over orders - starting from order 2 and ending at order N-1 */
    for (i=1; i<spec->nor-1; i++) {
      /* Read seeing for this order */
      sprintf(card,"HIERARCH ESO QC ORD%d OBJ FWHM",i+1);
      if (fits_read_key(infits,TDOUBLE,card,&(spec->or[i].seeing),comment,&status))
	errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.",card,spec->file);
      /* Convert from spatial pixels to arceseconds */
      spec->or[i].seeing*=spec->pixscal*(double)spec->biny;
    }
    /* If values first and last orders are not present then replace
       them with the values for order 2 and N-1 respectively */
    sprintf(card,"HIERARCH ESO QC ORD%d OBJ FWHM",1);
    if (fits_read_key(infits,TDOUBLE,card,&(spec->or[0].seeing),comment,&status)) {
      status=0; spec->or[0].seeing=spec->or[1].seeing;
    }
    sprintf(card,"HIERARCH ESO QC ORD%d OBJ FWHM",spec->nor);
    if (fits_read_key(infits,TDOUBLE,card,&(spec->or[spec->nor-1].seeing),comment,
		      &status)) {
      status=0; spec->or[spec->nor-1].seeing=spec->or[spec->nor-2].seeing;
    }
    /* Calculate the median seeing in this spectrum for writing to spectrum structure */
    /* Allocate memory for temporary seeing array */
    if ((seeing=darray(spec->nor))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for temp.\n\
\tseeing array of size %d for file\n\t%s",spec->nor,spec->file);
    for (i=0,j=0; i<spec->nor; i++) if (spec->or[i].seeing>0.0) seeing[j++]=spec->or[i].seeing;
    stat.med=0.0; if (j>0 && !median(seeing,NULL,j,&stat,0))
      errormsg("UVES_r2Dspec(): Error returned from median() while\n\
\tanalysing seeing data from file\n\t%s",spec->file);
    spec->seeing[2]=stat.med;
    /* Clean up */
    free(seeing);
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
    if (!naxis) errormsg("UVES_r2Dspec_ESOmer(): Couldn't find input image extension\n\
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
      errormsg("UVES_r2Dspec(): Cannot open FITS file\n\t%s",spec->thfile);
    /* For now, advance one HDU in CPL-written files */
    if (spec->fvers==1 && fits_movrel_hdu(infits,1,&hdutype,&status))
      errormsg("UVES_r2Dspec(): Cannot move to second HDU\n\
\tin FITS file\n\t%s.",spec->thfile);

    /* Read in archival filename */
    if (!spec->fvers &&
	fits_read_key(infits,TSTRING,"ARCFILE",spec->tharfile,comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","ARCFILE",spec->thfile);
    else if (spec->fvers==1) {
      if (fits_read_key(infits,TSTRING,"MJD-OBS",card,comment,&status))
	errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->thfile);
      if ((cptr=strchr(comment,'('))==NULL)
	errormsg("UVES_r2Dspec(): Cannot find start of\n\
\tUT time stamp in comment of header card\n\t%s in file\n\t%s","MJD-OBS",spec->thfile);
      sprintf(spec->tharfile,"UVES.%s",cptr+1);
      if ((cptr=strchr(spec->tharfile,')'))==NULL)
	errormsg("UVES_r2Dspec(): Cannot find end of\n\
\tUT time stamp in comment of header card\n\t%s in file\n\t%s","MJD-OBS",spec->thfile);
      sprintf(cptr,"%s",".fits");
    }

    /* Get modified julian day */
    if (fits_read_key(infits,TDOUBLE,"MJD-OBS",&(spec->wc_jd),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->thfile);
    /* Convert to Julian day */
    spec->wc_jd+=2400000.5;

    /* Find the temperature in each arm and the atmospheric pressure */
    if (fits_read_key(infits,TSTRING,"HIERARCH ESO INS PATH",card,comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO INS PATH",spec->thfile);
    if (strstr(card,"RED")!=NULL || strstr(card,"red")!=NULL)
      sprintf(card,"HIERARCH ESO INS TEMP2 MEAN");
    else sprintf(card,"HIERARCH ESO INS TEMP1 MEAN");
    if (fits_read_key(infits,TDOUBLE,card,&(spec->wc_temp),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.",card,spec->thfile);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SENS26 MEAN",
		      &(spec->wc_pres),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SENS26 MEAN",spec->thfile);

    /* For now we have to move back one HDU if dealing with CPL-written files */
    if (spec->fvers==1 && fits_movrel_hdu(infits,-1,&hdutype,&status))
      errormsg("UVES_r2Dspec(): Cannot move to second HDU\n\
\tin FITS file\n\t%s.",spec->thfile);

    /* Get image dimensions */
    if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
      errormsg("UVES_r2Dspec(): Couldn't get image dimensions for\n\
\tfile %s",spec->thfile);
    if (!naxis) errormsg("UVES_r2Dspec(): Couldn't find input image extension\n\
\tfor file %s",spec->thfile);

    /* Some weak checks that this really is the right ThAr array */
    if (naxes[0]!=npix)
      errormsg("UVES_r2Dspec(): Mismatch in array sizes between\n\
\tflux and ThAr files\n\t%s &\n\t%s",spec->file,spec->thfile);
    if (naxes[1]!=spec->nor)
      errormsg("UVES_r2Dspec(): Mismatch in number of orders between\n\
\tflux and ThAr files\n\t%s &\n\t%s",spec->file,spec->thfile);

    /* Read in ThAr information */
    for (i=0,first=1; i<spec->nor; i++,first+=npix) {
      if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
			spec->or[i].th,&anynul,&status))
	errormsg("UVES_r2Dspec(): Cannot read ThAr flux array for order\n\
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
	if (spec->or[i].fl[j]==0.0) spec->or[i].st[j]=RCLIP;
      /* Determine where the useful part of the spectrum begins and ends */
      j=0; while (j<spec->or[i].np && spec->or[i].st[j]<0) j++;
      k=j; j=spec->or[i].np-1; while (j>=0 && spec->or[i].st[j]<0) j--;
      spec->or[i].nuse=j-k+1;
    }
    if (spec->or[i].nuse<MINUSE) {
      /* Free the allocated memory for the raw arrays */
      if (!UVES_memspec(spec,par,-i-1,0))
	errormsg("UVES_r2Dspec(): Error returned from UVES_memspec() when\n\
\tattempting to free memory from order %d of file\n\t%s",i+1,spec->file);
      /*
      free(spec->or[i].vhwl); free(spec->or[i].vhrwl); free(spec->or[i].fl);
      free(spec->or[i].er); free(spec->or[i].res); free(spec->or[i].st);
      if (par->thar==1) {
	free(spec->or[i].th); free(spec->or[i].ter); free(spec->or[i].tst);
      }
      */
    }
  }

  /* Attempt to open wavelength solution fits file */
  if (fits_open_file(&infits,spec->wlfile,READONLY,&status))
      errormsg("UVES_r2Dspec(): Cannot open FITS file\n\t%s",spec->wlfile);
  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec(): Cannot find number of HDUs in file\n\
\t%s",spec->wlfile);
  if (!spec->fvers && hdunum!=2)
    errormsg("UVES_r2Dspec(): Number of HDUs is %d instead of %d\n\
\tin file %s",hdunum,2,spec->wlfile);
  else if (spec->fvers==1 && hdunum!=10)
    errormsg("UVES_r2Dspec(): Number of HDUs is %d instead of %d\n\
\tin file %s",hdunum,10,spec->wlfile);
  /* Read in archival filename */
  if (!spec->fvers &&
      fits_read_key(infits,TSTRING,"ARCFILE",spec->wlarfile,comment,&status))
    errormsg("UVES_r2Dspec(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","ARCFILE",spec->wlfile);
  else if (spec->fvers==1) {
    if (fits_read_key(infits,TSTRING,"MJD-OBS",card,comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->wlfile);
    if ((cptr=strchr(comment,'('))==NULL)
      errormsg("UVES_r2Dspec(): Cannot find start of\n\
\tUT time stamp in comment of header card\n\t%s in file\n\t%s","MJD-OBS",spec->wlfile);
    sprintf(spec->wlarfile,"UVES.%s",cptr+1);
    if ((cptr=strchr(spec->wlarfile,')'))==NULL)
      errormsg("UVES_r2Dspec(): Cannot find end of\n\
\tUT time stamp in comment of header card\n\t%s in file\n\t%s","MJD-OBS",spec->wlfile);
    sprintf(cptr,"%s",".fits");
  }
  /* Move to relevant HDU */
  if (!spec->fvers && fits_movrel_hdu(infits,1,&hdutype,&status))
    errormsg("UVES_r2Dspec(): Could not move to second HDU\n\
\tin file %s",spec->wlfile);
  else if (spec->fvers==1 && fits_movrel_hdu(infits,5,&hdutype,&status))
    errormsg("UVES_r2Dspec(): Could not move to sixth HDU\n\
\tin file %s",spec->wlfile);

  /* Read in wavelength solution and, for CPL files, ThAr header
     information and line parameters */
  /* Check HDU type */
  if (!spec->fvers && hdutype!=BINARY_TBL)
    errormsg("UVES_r2Dspec(): Extension %d not a binary table\n\
\tin file\n\t%s",2,spec->wlfile);
  else if (spec->fvers==1 && hdutype!=BINARY_TBL)
    errormsg("UVES_r2Dspec(): Extension %d not a binary table\n\
\tin file\n\t%s",6,spec->wlfile);
  if (!spec->fvers) {
    /* Search HISTORY records for wavelength information */
    *inclist=search;
    while (!fits_find_nextkey(infits,inclist,1,inclist,0,card,&status) &&
	   strncmp(card,"HISTORY  'COEFD'",16));
    if (strncmp(card,"HISTORY  'COEFD'",16))
      errormsg("UVES_r2Dspec(): Cannot find header record beginning\n\
\twith %s in file \n\t%s","HISTORY  'COEFD'",spec->wlfile);
    while ((cptr=strchr(card,','))!=NULL) *cptr=' ';
    if (sscanf(card,"%s %s %s %d %d %s",dummy,dummy,dummy,&idum,&ncoefd,dummy)!=6)
      errormsg("UVES_r2Dspec(): Cannot read number of wavelength\n\
\tcoefficients available from file \n\t%s",spec->wlfile);
    /* Assign number of wavelegth coefficients */
    for (i=1,spec->or[0].nwpol=ncoefd/spec->nor; i<spec->nor; i++)
      spec->or[i].nwpol=spec->or[0].nwpol;
    if (spec->or[0].nwpol>NWPOL)
      errormsg("UVES_r2Dspec(): Predicted number of wavelength\n\
\tcoefficients per order, %d, is greater than NWPOL(=%d) in file\n\t%s.\n\
\tIncrease NWPOL and recompile\n\t%s",spec->or[0].nwpol,NWPOL,spec->wlfile);
    if (ncoefd!=spec->nor*spec->or[0].nwpol)
      errormsg("UVES_r2Dspec(): Number of wavelength coefficients\n\
\tavailable is not equal to predicted number of coefficients\n\
\tper order, %d, times the number of orders, %d, in file\n\t%s",spec->or[0].nwpol,
	       spec->nor,spec->wlfile);
    /* Allocate memory for coefficient array */
    if ((coefd=darray(ncoefd))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for coefd\n\
\tarray of size %d for file\n\t%s",ncoefd,spec->wlfile);
    /* Read in all coefficients */
    for (i=0,j=0; i<ncoefd/3; i++,j+=3) {
      if (fits_find_nextkey(infits,inclist,1,inclist,0,card,&status))
	errormsg("UVES_r2Dspec(): Cannot find wavelength information\n\
\trecord %d in file\n\t%s",i+1,spec->wlfile);
      if (sscanf(card,"%s %lf %lf %lf",dummy,&(coefd[j]),&(coefd[j+1]),
		 &(coefd[j+2]))!=4)
	errormsg("UVES_r2Dspec(): Wrong format in wavelength\n\
information block, line %d, in file\n\t%s",i+1,spec->wlfile);
    }
    if (j<spec->nor) {
      if (fits_find_nextkey(infits,inclist,1,inclist,0,card,&status))
	errormsg("UVES_r2Dspec(): Cannot find wavelength information\n\
\trecord %d in file\n\t%s",i+2,spec->wlfile);
      if (sscanf(card,"%s %lf %lf",dummy,&(coefd[j]),&(coefd[j+1]))!=ncoefd-j+1)
	errormsg("UVES_r2Dspec(): Wrong format in wavelength\n\
information block, line %d, in file\n\t%s",i+2,spec->wlfile);
    }
    /* Check that table is 2-dimensional */
    if (fits_read_key(infits,TINT,"NAXIS",&naxis,comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS",spec->wlfile);
    if (naxis!=2) errormsg("UVES_r2Dspec(): The binary table in file\n\t%s\n\
\tis %d-dimensional. It should be 2-dimensional",spec->wlfile,naxis);
    /* Find number of axes to be read in */
    if (fits_get_num_cols(infits,&naxis,&status))
      errormsg("UVES_r2Dspec(): Cannot read number of columns\n\
\t%s in FITS file\n\t%s.","NAXIS2",spec->wlfile);
    if (naxis<17) errormsg("UVES_r2Dspec(): The binary table in file\n\t%s\n\
\thas %d columns. It should have at least 17",spec->wlfile,naxis);
    /* Find number of rows to be read in */
    if (fits_get_num_rows(infits,&nrows,&status))
      errormsg("UVES_r2Dspec(): Cannot read number of rows\n\
\t%s in FITS file\n\t%s.","NAXIS2",spec->wlfile);
    /* Allocate memory for resolution array */
    if ((resol=darray(nrows))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for resol\n\
\tarray of size %d for file\n\t%s",nrows,spec->wlfile);
    /* Read in resolution data */
    if (fits_read_col(infits,TDOUBLE,17,1,1,nrows,&nulval,resol,&anynul,&status))
      errormsg("UVES_r2Dspec(): Cannot read resolution array in file\n\t%s",
	       spec->wlfile);
    /* Allocate the coefficients to their respective orders */
    for (i=0; i<spec->nor; i++) {
      for (j=0,k=i*spec->or[0].nwpol; j<spec->or[0].nwpol; j++,k++)
	spec->or[i].wpol[j]=coefd[k];
    }
    /* Determine median resolution */
    i=0; while (i<nrows && resol[i]>0.0) i++;
    if (!median(resol,NULL,i,&stat,0))
      errormsg("UVES_r2Dspec(): Error returned from median() while\n\
\tanalysing resolution data from file\n\t%s",spec->wlfile);
    spec->arcfwhm=C_C_K/stat.med;
    /* Clean up */
    free(coefd); free(resol);
  } else if (spec->fvers==1) {
    /* Polynomial coefficients for each echelle order need to be
       constructed from 2D fit coefficients */
    /* Find number of axes to be read in */
    if (fits_get_num_cols(infits,&naxis,&status))
      errormsg("UVES_r2Dspec(): Cannot read number of columns\n\
\t%s in FITS file\n\t%s.","TFIELDS",spec->wlfile);
    if (naxis!=3) errormsg("UVES_r2Dspec(): The binary table in file\n\t%s\n\
\thas %d columns. It should have 3",spec->wlfile,naxis);
    /* Find number of rows to be read in */
    if (fits_get_num_rows(infits,&nrows,&status))
      errormsg("UVES_r2Dspec(): Cannot read number of rows\n\
\t%s in FITS file\n\t%s.","NAXIS2",spec->wlfile);
    /* Allocate memory to hold index and coefficient arrays */
    if ((coefd=darray(nrows))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for coefd\n\
\tarray of size %d",nrows);
    if ((ord1=iarray(nrows))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for ord1\n\
\tarray of size %d",nrows);
    if ((ord2=iarray(nrows))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for ord2\n\
\tarray of size %d",nrows);
    /* Find and read in the diffraction order numbers from HISTORY cards */
    *inclist=search;
    while (!fits_find_nextkey(infits,inclist,1,inclist,0,card,&status) &&
	   strncmp(card,"HISTORY FABSORD",15));
    if (strncmp(card,"HISTORY FABSORD",15))
      errormsg("UVES_r2Dspec(): Cannot find header record beginning\n\
\twith %s in file \n\t%s","HISTORY FABSORD",spec->wlfile);
    if (sscanf(card,"%s %s %d",dummy,dummy,&ms)!=3)
      errormsg("UVES_r2Dspec(): Cannot read starting diffraction\n\
\torder number from file \n\t%s",spec->wlfile);
    while (!fits_find_nextkey(infits,inclist,1,inclist,0,card,&status) &&
	   strncmp(card,"HISTORY LABSORD",15));
    if (strncmp(card,"HISTORY LABSORD",15))
      errormsg("UVES_r2Dspec(): Cannot find header record beginning\n\
\twith %s in file \n\t%s","HISTORY LABSORD",spec->wlfile);
    if (sscanf(card,"%s %s %d",dummy,dummy,&me)!=3)
      errormsg("UVES_r2Dspec(): Cannot read ending diffraction\n\
\torder number from file \n\t%s",spec->wlfile);
    /* Check that diffraction order numbers are consistent with number
       of axis read */
    if (abs(ms-me+1)!=spec->nor)
      errormsg("UVES_r2Dspec(): Number of orders read in is &d\n\
\tbut labelling of diffraction orders goes from %d to %d\n\
\t(=%d orders) in file\n\
\t%s.\n\tThis inconsistency must be resolved for even remotely reliable\n\
\twavelength solutions to be applied!",spec->nor,ms,me,abs(ms-me+1),spec->wlfile);
    /* Define direction of increasing diffraction order */
    mslope=1; if (me<ms) mslope=-1;
    /* Find and read in the index and coefficient arrays */
    if (fits_get_colnum(infits,CASEINSEN,"ORDER1",&col,&status))
      errormsg("UVES_r2Dspec(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","ORDER1",spec->wlfile);
    if (fits_read_col(infits,TINT,col,1,1,nrows,&nulval,ord1,&anynul,&status))
      errormsg("UVES_r2Dspec(): Cannot read column '%s'\n\
\tin binary table in FITS file\n\t%s","ORDER1",spec->wlfile);
    if (fits_get_colnum(infits,CASEINSEN,"ORDER2",&col,&status))
      errormsg("UVES_r2Dspec(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","ORDER2",spec->wlfile);
    if (fits_read_col(infits,TINT,col,1,1,nrows,&nulval,ord2,&anynul,&status))
      errormsg("UVES_r2Dspec(): Cannot read column '%s'\n\
\tin binary table in FITS file\n\t%s","ORDER2",spec->wlfile);
    if (fits_get_colnum(infits,CASEINSEN,"COEFF",&col,&status))
      errormsg("UVES_r2Dspec(): Cannot find column named '%s'\n\
\tin binary table in FITS file\n\t%s","COEFF",spec->wlfile);
    if (fits_read_col(infits,TDOUBLE,col,1,1,nrows,&nulval,coefd,&anynul,&status))
      errormsg("UVES_r2Dspec(): Cannot read column '%s'\n\
\tin binary table in FITS file\n\t%s","COEFF",spec->wlfile);
    /* Go through coefficient array to get offsets and check that
       normalizations haven't been used */
    for (i=0; i<6; i++) {
      /* Check that normalizations have not been used */
      if (i<2) spec->or[0].wpol[NWPOL-2+i]=coefd[i];
      else if (i==2) mnorm=coefd[i];
      else if (i>2 && coefd[i]!=1.0)
	errormsg("UVES_r2Dspec(): Non-unity normalization parameters\n\
\tin lines %d to %d of 'COEFF' column in binary table\n\
\tin FITS file\n\t%s",4,6,spec->wlfile);
    }
    /* Determine the degrees of the polynomials */
    /* Current algorithm is to go through array in reverse and find
       the first time the two index arrays are equal. This assumes
       that the degree in the spectral and spatial directions is the
       same */ 
    for (i=nrows-1; i>=6; i--) if (ord1[i]==ord2[i]) break;
    spec->or[0].nwpol=No=ord1[i]+1;
    /* Check that we have enough space for coefficients and fill
       number and offsets for other orders */
    if (spec->or[0].nwpol>NWPOL-2)
      errormsg("UVES_r2Dspec(): Order of wavelength calibration\n\
\tpolynomial in spectral direction (=%d) must be <= NWPOL-2=%d for\n\
\tUVES CPL-reduced data. Try increasing NWPOL to %d and recompiling",spec->or[0].nwpol,
	     NWPOL-2,spec->or[0].nwpol+2);
    for (i=1; i<spec->nor; i++) {
      spec->or[i].nwpol=spec->or[0].nwpol;
      spec->or[i].wpol[NWPOL-2]=spec->or[0].wpol[NWPOL-2];
      spec->or[i].wpol[NWPOL-1]=spec->or[0].wpol[NWPOL-1];
    }
    /* Loop to calculate coefficients for each order */
    for (i=0,m=ms; i<spec->nor; i++,m+=mslope) {
      spec->or[i].wpol[NWPOL-2]/=(double)m;
      dm=(double)m-mnorm;
      for (j=0; j<spec->or[i].nwpol; j++) {
	for (k=6,l=0,spec->or[i].wpol[j]=0.0; k<nrows; k++) {
	  if (ord1[k]==j) {
	    spec->or[i].wpol[j]+=coefd[k]*pow(dm,(double)ord2[k]); l++;
	  }
	  if (l==No) break;
	}
	spec->or[i].wpol[j]/=(double)m;
      }
    }
    /* Clean up */
    free(coefd); free(ord1); free(ord2);
    /* Move back to first HDU */
    if (fits_movabs_hdu(infits,1,&hdutype,&status))
      errormsg("UVES_r2Dspec(): Could not move to first HDU in file\n\t%s",
	       spec->wlfile);    
    /** Read in relevant header cards for CPL-reduced wavelength calibrations */
    /* Get modified julian day */
    if (fits_read_key(infits,TDOUBLE,"MJD-OBS",&(spec->wc_jd),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD-OBS",spec->wlfile);
    /* Convert to Julian day */
    spec->wc_jd+=2400000.5;
    /* Get date of observation and convert to year, month and day */
    if (fits_read_key(infits,TSTRING,"DATE-OBS",date,comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DATE-OBS",spec->wlfile);
    if (strlen(date)>10) {
      /* This means that the UT is appended to the date */
      if (sscanf(date,"%d-%d-%dT%lf:%lf:%lf",&(spec->wc_year),&(spec->wc_month),
		 &(spec->wc_day),&hour,&min,&sec)!=6)
	errormsg("UVES_r2Dspec(): Cannot read format of keyword\n\
%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->wlfile);
      /* Convert hours, mins and secs into the UT */
      spec->wc_ut=sec/3600.0+min/60.0+hour;
    } else {
      /* Only the date should be given */
      if (sscanf(date,"%d-%d-%d",&(spec->wc_year),&(spec->wc_month),&(spec->wc_day))!=3)
	errormsg("UVES_r2Dspec(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->wlfile);
      /* Get universal time from TM-START or UTC */ 
      if (fits_read_key(infits,TDOUBLE,"TM-START",&(spec->wc_ut),comment,&status)) {
	if (fits_read_key(infits,TDOUBLE,"UTC",&(spec->wc_ut),comment,&status))
	  errormsg("UVES_r2Dspec(): Cannot read values of header cards\n\
\t%s & %s from FITS file\n\t%s.","TM-START","UTC",spec->wlfile);
      }
      /* Convert to hours */
      spec->wc_ut/=3600.0;
    }
    /* Get ThAr exposure time */
    if (fits_read_key(infits,TDOUBLE,"EXPTIME",&(spec->wc_etime),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","EXPTIME",spec->wlfile);
    /* Read in ThAr slit width and encoder step (can, in principle, be different to
       that used for object) */
    if (spec->cwl<UVESBORR) {
      if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SLIT2 WID",&(spec->wc_sw),
			comment,&status))
	errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SLIT2 WID",spec->wlfile);
      if (fits_read_key(infits,TINT,"HIERARCH ESO INS GRAT1 ENC",&(spec->wc_encoder),
			comment,&status)) {
	errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS GRAT1 ENC",spec->wlfile);
      }
    } else {
      if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SLIT3 WID",&(spec->wc_sw),
			comment,&status))
	errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SLIT3 WID",spec->wlfile);
      if (fits_read_key(infits,TINT,"HIERARCH ESO INS GRAT2 ENC",&(spec->wc_encoder),
			comment,&status)) {
	errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS GRAT2 ENC",spec->wlfile);
      }
    }
    /* Find the temperature in each arm and the atmospheric pressure */
    if (fits_read_key(infits,TSTRING,"HIERARCH ESO INS PATH",card,comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ESO INS PATH",spec->wlfile);
    if (strstr(card,"RED")!=NULL || strstr(card,"red")!=NULL)
      sprintf(card,"HIERARCH ESO INS TEMP2 MEAN");
    else sprintf(card,"HIERARCH ESO INS TEMP1 MEAN");
    if (fits_read_key(infits,TDOUBLE,card,&(spec->wc_temp),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.",card,spec->wlfile);
    if (fits_read_key(infits,TDOUBLE,"HIERARCH ESO INS SENS26 MEAN",
		      &(spec->wc_pres),comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card\n\
\t%s from FITS file %s.","HIERARCH ESO INS SENS26 MEAN",spec->wlfile);
    /** Read in ThAr line set information **/
    /* Read in CCD binning values */
    if (fits_read_key(infits,TINT,"HIERARCH ESO DET WIN1 BINY",&(spec->ts.binx),
		      comment,&status))
      errormsg("UVES_r2Dspec(): Cannot read value of header card %s from FITS file\n\
\t%s.","HIERARCH ESO DET WIN1 BINY",spec->wlfile);
    /* Move to HDU containing ThAr information */
    if (fits_movrel_hdu(infits,4,&hdutype,&status))
      errormsg("UVES_r2Dspec(): Cannot move to extension %d (HDU %d) in FITS file\n\
\t%s.",4,5,spec->wlfile);
    /* Check HDU type */
    if (hdutype!=BINARY_TBL)
      errormsg("UVES_r2Dspec(): Extension %d (HDU %d) not a binary table in file\n\
\t%s",4,5,spec->wlfile);
    /* Find number of axes to be read in */
    if (fits_get_num_cols(infits,&naxis,&status))
      errormsg("UVES_r2Dspec(): Cannot read number of columns %s in FITS file\n\
\t%s","TFIELDS",spec->wlfile);
    if (naxis!=26) errormsg("UVES_r2Dspec(): The binary table in file\n\t%s\n\
\thas %d columns. It should have %d",spec->wlfile,naxis,26);
    /* Find number of rows to be read in */
    if (fits_get_num_rows(infits,&nrows,&status))
      errormsg("UVES_r2Dspec(): Cannot read number of rows %s in FITS file\n\t%s",
	       "NAXIS2",spec->wlfile);
    spec->ts.n=nrows;
    /* Allocate memory for ThAr line set */
    if ((spec->ts.x=darray(spec->ts.n))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for spec->ts.x\n\
\tarray of length %d",spec->ts.n);
    if ((spec->ts.w=darray(spec->ts.n))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for spec->ts.w\n\
\tarray of length %d",spec->ts.n);
    if ((spec->ts.dis=darray(spec->ts.n))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for spec->ts.dis\n\
\tarray of length %d",spec->ts.n);
    if ((spec->ts.wlf=darray(spec->ts.n))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for spec->ts.wlf\n\
\tarray of length %d",spec->ts.n);
    if ((spec->ts.resid=darray(spec->ts.n))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for spec->ts.resid\n\
\tarray of length %d",spec->ts.n);
    if ((spec->ts.wlsf=darray(spec->ts.n))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for spec->ts.wlsf\n\
\tarray of length %d",spec->ts.n);
    if ((spec->ts.stp=iarray(spec->ts.n))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for spec->ts.stp\n\
\tarray of length %d",spec->ts.n);
    /* Find and read in relevant columns */
    /* X Position in pixels */
    if (fits_get_colnum(infits,CASEINSEN,"X",&col,&status))
      errormsg("UVES_r2Dspec(): Cannot find column '%s' in binary table in file\n\t%s",
	       "X",spec->wlfile);
    if (fits_read_col(infits,TDOUBLE,col,1,1,nrows,&nulval,spec->ts.x,&anynul,&status))
      errormsg("UVES_r2Dspec(): Cannot read column '%s' in binary table in file\n\t%s",
	       "X",spec->wlfile);
    /* Width in pixels */
    if (fits_get_colnum(infits,CASEINSEN,"Xwidth",&col,&status))
      errormsg("UVES_r2Dspec(): Cannot find column '%s' in binary table in file\n\t%s",
	       "Xwidth",spec->wlfile);
    if (fits_read_col(infits,TDOUBLE,col,1,1,nrows,&nulval,spec->ts.w,&anynul,&status))
      errormsg("UVES_r2Dspec(): Cannot read column '%s' in binary table in file\n\t%s",
	       "Xwidth",spec->wlfile);
    /* Convert to FWHM */
    for (i=0; i<spec->ts.n; i++) spec->ts.w[i]*=C_FWHMSIG;
    /* Residual in Ang. */
    if (fits_get_colnum(infits,CASEINSEN,"Residual",&col,&status))
      errormsg("UVES_r2Dspec(): Cannot find column '%s' in binary table in file\n\t%s",
	       "Residual",spec->wlfile);
    if (fits_read_col(infits,TDOUBLE,col,1,1,nrows,&nulval,spec->ts.resid,&anynul,&status))
      errormsg("UVES_r2Dspec(): Cannot read column '%s' in binary table in file\n\t%s",
	       "Residual",spec->wlfile);
    /* Dispersion in Ang./pix */
    if (fits_get_colnum(infits,CASEINSEN,"Pixel",&col,&status))
      errormsg("UVES_r2Dspec(): Cannot find column '%s' in binary table in file\n\t%s",
	       "Pixel",spec->wlfile);
    if (fits_read_col(infits,TDOUBLE,col,1,1,nrows,&nulval,spec->ts.dis,&anynul,
		      &status))
      errormsg("UVES_r2Dspec(): Cannot read column '%s' in binary table in file\n\t%s",
	       "Pixel",spec->wlfile);
    /* Fitted wavelength of line [Ang.] */
    if (fits_get_colnum(infits,CASEINSEN,"WaveC",&col,&status))
      errormsg("UVES_r2Dspec(): Cannot find column '%s' in binary table in file\n\t%s",
	       "WaveC",spec->wlfile);
    if (fits_read_col(infits,TDOUBLE,col,1,1,nrows,&nulval,spec->ts.wlf,&anynul,
		      &status))
      errormsg("UVES_r2Dspec(): Cannot read column '%s' in binary table in file\n\t%s",
	       "WaveC",spec->wlfile);
    /* Whether line was used in polynomial solution */
    if (fits_get_colnum(infits,CASEINSEN,"NLinSol",&col,&status))
      errormsg("UVES_r2Dspec(): Cannot find column '%s' in binary table in file\n\t%s",
	       "NLinSol",spec->wlfile);
    if (fits_read_col(infits,TINT,col,1,1,nrows,&nulval,spec->ts.stp,&anynul,&status))
      errormsg("UVES_r2Dspec(): Cannot read column '%s' in binary table in file\n\t%s",
	       "NLinSol",spec->wlfile);
    /* Find median FWHM resolution (in km/s) by constructing a temporary array */
    /* Allocate memory for resolution array */
    if ((resol=darray(spec->ts.n))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for resol\n\
\tarray of size %d for file\n\t%s",spec->ts.n,spec->wlfile);
    for (i=0; i<spec->ts.n; i++) resol[i]=spec->ts.w[i]*spec->ts.dis[i]/spec->ts.wlf[i];
    stat.med=0.0; if (spec->ts.n>0 && !median(resol,spec->ts.stp,spec->ts.n,&stat,0))
      errormsg("UVES_r2Dspec(): Error returned from median() while\n\
\tanalysing resolution data from file\n\t%s",spec->wlfile);
    spec->arcfwhm=C_C_K*stat.med;
    /* Clean up */
    free(resol);
    /* Find RMS residual (in m/s) by constructing a temporary array */
    /* Allocate memory for residual array */
    if ((resid=darray(spec->ts.n))==NULL)
      errormsg("UVES_r2Dspec(): Cannot allocate memory for resid\n\
\tarray of size %d for file\n\t%s",spec->ts.n,spec->wlfile);
    for (i=0; i<spec->ts.n; i++) resid[i]=spec->ts.resid[i]/spec->ts.wlf[i];
    stat.rms=0.0;
    if (spec->ts.n>0 && !stats(resid,NULL,NULL,NULL,spec->ts.stp,spec->ts.n,0,&stat))
      errormsg("UVES_r2Dspec(): Error returned from stats() while\n\
\tcalculating RMS of residuals from file\n\t%s",spec->wlfile);
    spec->wc_resid=C_C*stat.rms;
    /* Clean up */
    free(resid);
    /* Count lines actually used in the polynomial solution */
    for (i=0,spec->ts.np=0; i<spec->ts.n; i++) if (spec->ts.stp[i]==1) spec->ts.np++;
    /* Use ThAr line parameters and seeing information to model the
       resolution as a function of position along each order */
    if (!UVES_model_resol(spec)) {
      nferrormsg("UVES_r2Dspec(): Error returned from UVES_model_resol()\n\
\twhen computing resolution arrays for echelle orders in file\n\t%s",spec->file);
      return 0;
    }
  }
  /* Close wavelength solution fits file */
  fits_close_file(infits,&status);

  /* Make sure this spectrum will be combined into combined spectrum later */
  spec->comb=1;

  return 1;

}
