/****************************************************************************
* Read in a 1-dimensional spectrum, either from a FITS file or from an
* ASCII file
****************************************************************************/

#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <longnam.h>
#include "UVES_popler.h"
#include "const.h"
#include "file.h"
#include "memory.h"
#include "error.h"

int UVES_r1Dspec(cspectrum *cspec, params *par) {

  double   version=0.0,crval=0.0,cdelt=0.0,hdisp=0.0;
  double   nulval=0.0;
  double   ddum=0.0;
  long     nrows=0;
  int      hdutype=0,hdunum=0,status=0,first=1,naxes=0,anynul=0;
  int      idx=0,crpix=0,dcflag=0,ESO_PHASE3=0,UVES_SQUAD_PHASE3=0;
  int      col=0;
  int      i=0;
  char     buffer[LNGSTRLEN]="\0";
  char     comment[FLEN_COMMENT]="\0";
  char     dummy[FLEN_KEYWORD]="\0";
  char     search[FLEN_CARD]="HISTORY\0",card[FLEN_CARD]="\0";
  char     *inclist[1];
  FILE     *data_file=NULL;
  fitsfile *infits;

  /** Open input file, see if it's a FITS file and branch as a result **/
  if ((data_file=faskropen("FITS or ASCII file name?",cspec->file,5))==NULL)
    errormsg("UVES_r1Dspec(): Can not open file %s",cspec->file);
  /* Read in first line from input file to decide its nature */
  idx++; if (fgets(buffer,LNGSTRLEN,data_file)==NULL) {
    fclose(data_file);
    errormsg("UVES_r1Dspec(): Problem reading file %s on line %d",cspec->file,idx);
  }
  if (!strncmp(buffer,"SIMPLE  =",8)) {
    /* Assume that this means that this is a FITS file */
    /* Close file */
    fclose(data_file);
    /* Re-open input file as FITS file */
    if (fits_open_file(&infits,cspec->file,READONLY,&status))
      errormsg("UVES_r1Dspec(): Cannot open FITS file\n\t%s",cspec->file);
    /* Check HDU type */
    if (fits_get_hdu_type(infits,&hdutype,&status))
      errormsg("UVES_r1Dspec(): Cannot get HDU type for file\n\t%s",cspec->file);
    if (hdutype!=IMAGE_HDU)
      errormsg("UVES_r1Dspec(): File not a FITS image: %s",cspec->file);
    /* Check number of HDUs */
    if (fits_get_num_hdus(infits,&hdunum,&status))
      errormsg("UVES_r1Dspec(): Cannot find number of HDUs in file\n\
\t%s",cspec->file);
    if (hdunum<1 || hdunum>7)
      errormsg("UVES_r1Dspec(): Number of HDUs is %d instead of a\n\
\tmaximum of %d in file\n\t%s",hdunum,7,cspec->file);
    /** See if the FITS file is one written out by UVES_popler **/
    /* Find the version number from the second line of the HISTORY statements */
    *inclist=search;
    while (!fits_find_nextkey(infits,inclist,1,inclist,0,card,&status) &&
	   strncmp(card,"HISTORY UVES_popler: Version",28));
    status=0;
    if (!strncmp(card,"HISTORY UVES_popler: Version",28)) {
      /* This seems to be a UVES_popler FITS file */
      /* Get version number */
      if (sscanf(card,"%s %s %s %lf*",dummy,dummy,dummy,&version)!=4)
	errormsg("UVES_r1Dspec(): I suspect that the file\n\t%s\n\
\twas written by UVES_popler but I cannot read the version number\n\
\tfrom keyword beginning with 'HISTORY UVES_popler: Version'",cspec->file);
      /* Find number of axes to be read in */
      if (fits_read_key(infits,TINT,"NAXIS2",&naxes,comment,&status))
	errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS2",cspec->file);
      if (naxes!=9) errormsg("UVES_r1Dspec(): Number of arrays to read in\n\
\tshould be %d but is NAXIS2=%d in file\n\t%s",9,naxes,cspec->file);
      /* Find number of pixels to read in */
      if (fits_read_key(infits,TINT,"NAXIS1",&(cspec->np),comment,&status))
	errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS1",cspec->file);
      /* Find index of starting pixel */
      if (fits_read_key(infits,TINT,"CRPIX1",&crpix,comment,&status))
	errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","CRPIX1",cspec->file);
      /* Find flag which determines whether spectrum is linear or log-linear */
      if (fits_read_key(infits,TINT,"DC-FLAG",&dcflag,comment,&status))
	errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","DC-FLAG",cspec->file);
      /* Find starting wavelength */
      if (fits_read_key(infits,TDOUBLE,"CRVAL1",&crval,comment,&status))
	errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","CRVAL1",cspec->file);
      if (crval<5.0) crval=pow(10.0,crval);
      /* Find wavelength spacing */
      if (fits_read_key(infits,TDOUBLE,"CDELT1",&cdelt,comment,&status)) {
	status=0;
	if (fits_read_key(infits,TDOUBLE,"CD1_1",&cdelt,comment,&status))
	  errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s or %s from FITS file\n\t%s.","CDELT1","CD1_1",cspec->file);
      }
      /* Set combined spectrum parameters depending on whether
	 spectrum is linear or log-linear */
      if (!dcflag) {
	par->linear=1;
	par->disp=cspec->dwl=cdelt; cspec->dv=0.0; cspec->flwl=crval-0.5*cdelt;
      }
      else {
	par->linear=0;
	cspec->dv=cdelt; par->disp=C_C_K*(pow(10.0,cspec->dv)-1.0); cspec->dwl=0.0;
	cspec->flwl=crval*pow(10.0,-0.5*cspec->dv);
      }
      /* Initialize the combined spectrum and set its wavelength scale */
      if (!UVES_init_cspec(cspec,par,0)) {
	fits_close_file(infits,&status);
	nferrormsg("UVES_r1Dspec(): Error returned from\n\tUVES_init_cspec()");
	return 0;
      }
      /* Read in normalized flux data */
      if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->no,
			&anynul,&status))
	errormsg("UVES_r1Dspec(): Cannot read normalized flux array in file\n\
\t%s",cspec->file);
      /* Read in normalzed error data */
      first+=cspec->np;
      if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->ne,
			&anynul,&status))
	errormsg("UVES_r1Dspec(): Cannot read normalized error array in file\n\
\t%s",cspec->file);
      /* Read in normalized expected fluctuation data */
      first+=cspec->np;
      if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->nf,
			&anynul,&status))
	errormsg("UVES_r1Dspec(): Cannot read normalized expected fluctuation\n\
\t array in file\n\t%s",cspec->file);
      /* Read in continuum data */
      first+=cspec->np;
      if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->co,
			&anynul,&status))
	errormsg("UVES_r1Dspec(): Cannot read continuum array in file\n\
\t%s",cspec->file);
      if (par->nocont==CNTNON) par->nocont=CNTSRD;
      /* Read in status data */
      first+=cspec->np;
      if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->fl,
			&anynul,&status))
	errormsg("UVES_r1Dspec(): Cannot read continuum array in file\n\
\t%s",cspec->file);
      /* Read in number of pixels before clip data */
      first+=cspec->np;
      if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->er,
			&anynul,&status))
	errormsg("UVES_r1Dspec(): Cannot read array of number of pixels\n\
\tbefore sigma clipping in file\n\t%s",cspec->file);
      /* Read in number of pixels after clip data */
      first+=cspec->np;
      if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->ef,
			&anynul,&status))
	errormsg("UVES_r1Dspec(): Cannot read array of number of pixels\n\
\tafter sigma clipping in file\n\t%s",cspec->file);
      /* Read in chisq. before clip data */
      first+=cspec->np;
      if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->csq,
			&anynul,&status))
	errormsg("UVES_r1Dspec(): Cannot read array of chisq.\n\
\tbefore sigma clipping in file\n\t%s",cspec->file);
      /* Read in chisq. after clip data */
      first+=cspec->np;
      if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->ccsq,
			&anynul,&status))
	errormsg("UVES_r1Dspec(): Cannot read array of chisq.\n\
\tafter sigma clipping in file\n\t%s",cspec->file);
      /* Set status array and un-normalized flux, error and expected
	 fluctuation arrays */
      for (i=0; i<cspec->np; i++) {
	cspec->st[i]=(NINT(cspec->fl[i])); cspec->ncb[i]=(NINT(cspec->er[i]));
	cspec->nccb[i]=(NINT(cspec->ef[i]));
	cspec->fl[i]=cspec->no[i]*cspec->co[i];
	cspec->er[i]=cspec->ne[i]*cspec->co[i];
	cspec->ef[i]=cspec->nf[i]*cspec->co[i];
      }
      /* Make sure there's a reasonable amount of data */
      if (cspec->np<MINUSE)
	errormsg("UVES_r1Dspec(): Less than MINUSE=%d data points in file\n\
\t%s.\n\tTry reducing MINUSE or zero-pad data file",MINUSE,cspec->file);
    } else {
      /* Not a UVES_popler FITS file so treat as a general spectral FITS file */
      /* First check to see whether this is a normal "multispec"
	 formatted file or if it's a tabular FITS file. Check is done
	 by simply seeing if number of pixels to be read in from
	 primary array is zero or not */
      /* Find number of axes to read in */
      if (fits_read_key(infits,TINT,"NAXIS",&naxes,comment,&status))
	errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS",cspec->file);
      if (naxes) {
	/* Check whether it's an SDSS file */
	if (fits_read_key(infits,TSTRING,"TELESCOP",card,comment,&status))
	  version=status=0;
	else version=(strstr(card,"SDSS 2.5-M")!=NULL) ? 1 : 0;
	/* This should be a normal "multispec" file */
	/* For the moment, just assume that if the file contains 3
	   columns they must be normalized flux, normalized error and
	   continuum, while if there's 2 or more than 3, just read
	   first two as flux and error */
	/* Find number of pixels to read in */
	if (fits_read_key(infits,TINT,"NAXIS1",&(cspec->np),comment,&status))
	  errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS1",cspec->file);
	/* Find number of axes to be read in */
	if (fits_read_key(infits,TINT,"NAXIS2",&naxes,comment,&status))
	  errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS2",cspec->file);
	warnmsg("UVES_r1Dspec(): Input file\n\t%s\n\
\twas not written by UVES_popler. Attempting to read file anyway.",cspec->file);
	if (naxes<2) errormsg("UVES_r1Dspec(): Only %d arrays in file\n\
\t%s.\n\tMust be at least 2 arrays, flux and error, to make sense",naxes,cspec->file);
	else if (naxes==2 || naxes>3)
	  warnmsg("UVES_r1Dspec(): There are %d data arrays to be read from file\n\
\t%s.\n\tAssuming that first two are flux and error, respectively",naxes,cspec->file);
	else
	  warnmsg("UVES_r1Dspec(): There are %d data arrays to be read from file\n\
\t%s.\n\tAssuming that the first 3 are normalized flux, normalized error\n \
\tand continuum, respectively",naxes,cspec->file);
	/* Find index of starting pixel */
	if (fits_read_key(infits,TINT,"CRPIX1",&crpix,comment,&status))
	  errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","CRPIX1",cspec->file);
	/* Find flag which determines whether spectrum is linear or log-linear */
	if (fits_read_key(infits,TINT,"DC-FLAG",&dcflag,comment,&status))
	  errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","DC-FLAG",cspec->file);
	/* Find starting wavelength */
	if (fits_read_key(infits,TDOUBLE,"CRVAL1",&crval,comment,&status))
	  errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","CRVAL1",cspec->file);
	if (crval<5.0) crval=pow(10.0,crval);
	/* Find wavelength spacing */
	if (fits_read_key(infits,TDOUBLE,"CDELT1",&cdelt,comment,&status)) {
	  status=0;
	  if (fits_read_key(infits,TDOUBLE,"CD1_1",&cdelt,comment,&status))
	    errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s or %s from FITS file\n\t%s.","CDELT1","CD1_1",cspec->file);
	}
	/* Set combined spectrum parameters depending on whether
	   spectrum is linear or log-linear */
	if (!dcflag) {
	  par->disp=cspec->dwl=cdelt; cspec->dv=0.0; cspec->flwl=crval-0.5*cdelt;
	}
	else {
	  cspec->dv=cdelt; par->disp=C_C_K*(pow(10.0,cspec->dv)-1.0); cspec->dwl=0.0;
	  cspec->flwl=crval*pow(10.0,-0.5*cspec->dv);
	}
	/* Initialize the combined spectrum and set its wavelength scale */
	if (!UVES_init_cspec(cspec,par,0)) {
	  fits_close_file(infits,&status);
	  nferrormsg("UVES_r1Dspec(): Error returned from\n\tUVES_init_cspec()");
	  return 0;
	}
	/* Read in flux data */
	if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->fl,
			  &anynul,&status))
	  errormsg("UVES_r1Dspec(): Cannot read flux array in file\n\
\t%s",cspec->file);
	/* Read in error data */
	first+=(version==1) ? 2*cspec->np : cspec->np;
	if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->er,
			  &anynul,&status))
	  errormsg("UVES_r1Dspec(): Cannot read error array in file\n\
\t%s",cspec->file);
	if (naxes==3) {
	  /* Read in continuum data */
	  first+=cspec->np;
	  if (fits_read_img(infits,TDOUBLE,first,cspec->np,&nulval,cspec->co,
			    &anynul,&status))
	    errormsg("UVES_r1Dspec(): Cannot read continuum array in file\n\
\t%s",cspec->file);
	  if (par->nocont==CNTNON) par->nocont=CNTSRD;
	  for (i=0; i<cspec->np; i++) {
	    cspec->no[i]=cspec->fl[i]; cspec->fl[i]*=cspec->co[i];
	    cspec->ne[i]=cspec->er[i]; cspec->er[i]*=cspec->co[i];
	  }
	} else {
	  for (i=0; i<cspec->np; i++) {
	    cspec->no[i]=cspec->fl[i]; cspec->ne[i]=cspec->er[i]; cspec->co[i]=1.0;
	  }
	}
	/* Initialize other arrays and check for bad pixels */
	for (i=0; i<cspec->np; i++) {
	  if (cspec->er[i]<=0.0) {
	    cspec->st[i]=RCLIP; cspec->no[i]=1.0; cspec->ncb[i]=cspec->nccb[i]=0;
	    cspec->ef[i]=cspec->er[i]=cspec->nf[i]=cspec->ne[i]=-INFIN;
	  } else {
	    cspec->st[i]=cspec->ncb[i]=cspec->nccb[i]=1;
	    cspec->ef[i]=cspec->er[i]; cspec->nf[i]=cspec->ne[i];
	  }
	  cspec->csq[i]=cspec->ccsq[i]=0.0;
	}
	/* Make sure there's a reasonable amount of data */
	if (cspec->np<MINUSE)
	  errormsg("UVES_r1Dspec(): Less than MINUSE=%d data points in file\n\
\t%s.\n\tTry reducing MINUSE or zero-pad data file",MINUSE,cspec->file);
      } else {
	/** This might be a tabular FITS file **/
	/* Check if this is a Phase 3 UVES_SQUAD file */
	if (fits_read_key(infits,TSTRING,"REFERENC",dummy,comment,&status)) status=0;
	else if (!strncmp(comment,"Reference for UVES_SQUAD DR",27)) UVES_SQUAD_PHASE3=1;
	/* Move to next HDU */
	if (fits_movrel_hdu(infits,1,&hdutype,&status))
	  errormsg("UVES_r1Dspec(): Cannot move to second HDU\n\tin FITS file\n\t%s.",
		   cspec->file);
	/* Check HDU type */
	if (hdutype!=BINARY_TBL)
	  errormsg("UVES_r2Dspec(): No primary array found in file\n\t%s\n\
\tso I'm looking for a binary table in the second FITS HDU.\n\
\tHowever, the second HDU is not a binary table",cspec->file);
	/* Check that table is 2-dimensional */
	if (fits_read_key(infits,TINT,"NAXIS",&naxes,comment,&status))
	  errormsg("UVES_r1Dspec(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS",cspec->file);
	if (naxes!=2) errormsg("UVES_r1Dspec(): The binary table in file\n\t%s\n\
\tis %d-dimensional. It should be 2-dimensional",cspec->file,naxes);
	/* Check if this is an ESO Phase 3 file */
	if (fits_read_key(infits,TSTRING,"VOPUB",card,comment,&status)) status=0;
	else if (!strncmp(card,"ESO/SAF",7)) ESO_PHASE3=1;
	/* Find number of axes to be read in */
        if (fits_get_num_cols(infits,&naxes,&status))
	  errormsg("UVES_r1Dspec(): Cannot read number of columns\n\
\tin FITS file\n\t%s.",cspec->file);
	/* If it's an ESO Phase 3 spectrum then check some basics */
	if (ESO_PHASE3) {
	  /* Check to make sure the basic set of columns can be read in later */
	  if (fits_get_colname(infits,CASEINSEN,"WAVE",card,&col,&status))
	    errormsg("UVES_r1Dspec(): Cannot find table column named %s\n\
\tin binary table in FITS file\n\t%s","WAVE_EL",cspec->file);
	  if (fits_get_colname(infits,CASEINSEN,"FLUX",card,&col,&status)) {
	    status=0;
	    if (fits_get_colname(infits,CASEINSEN,"FLUX_EL",card,&col,&status))
	      errormsg("UVES_r1Dspec(): Cannot find table column named %s\n\
\tin binary table in FITS file\n\t%s","FLUX or FLUX_EL",cspec->file);
	  }
	  if (fits_get_colname(infits,CASEINSEN,"ERR",card,&col,&status)) {
	    status=0;
	    if (fits_get_colname(infits,CASEINSEN,"FLUX_EL",card,&col,&status))
	      errormsg("UVES_r1Dspec(): Cannot find table column named %s\n\
\tin binary table in FITS file\n\t%s","ERR or ERR_EL",cspec->file);
	  }
	}
	/* If it's not an ESO Phase 3 spectrum then you have to assume some things */
	else if (naxes<3) errormsg("UVES_r1Dspec(): Only %d arrays in file\n\
\t%s.\n\tMust be at least 3 arrays, flux and error, to make sense",naxes,
			      cspec->file);
	else if (naxes==3)
	  warnmsg("UVES_r1Dspec(): There are %d data arrays to be read from file\n\
\t%s.\n\tAssuming that they are wavelength, flux and error, respectively",naxes,
		  cspec->file);
	else if (naxes>=4)
	  warnmsg("UVES_r1Dspec(): There are %d data arrays to be read from file\n\
\t%s.\n\tAssuming that the first is wavelength and the last 3 are normalized\n\
\tflux, normalized error and continuum, respectively",naxes,cspec->file);
	/* Find number of pixels to read in */
        if (fits_get_num_rows(infits,&nrows,&status))
	  errormsg("UVES_r1Dspec(): Cannot read number of rows\n\
\tin FITS file\n\t%s.",cspec->file);
	/* If there's a single row then seek the "NELEM" keyword for
	   clarification about the length of the arrays */
	if (nrows==1 || ESO_PHASE3) {
	  if (fits_read_key(infits,TLONG,"NELEM",&nrows,comment,&status))
	    errormsg("UVES_r1Dspec(): Number of rows is 1 and/or this is\n\
\tan ESO_PHASE3 file but the NELEM keyword cannot be found to clarify\n\
\tlength of data arrays in FITS file\n\t%s.",cspec->file);
	}
	cspec->np=(int)nrows;
	/* Initialize the combined spectrum but do not set its
	   wavelength scale */
	if (!UVES_init_cspec(cspec,par,1)) {
	  fits_close_file(infits,&status);
	  nferrormsg("UVES_r1Dspec(): Error returned from\n\tUVES_init_cspec()");
	  return 0;
	}
	/* Read in wavelength data */
	if (ESO_PHASE3) fits_get_colname(infits,CASEINSEN,"WAVE",card,&col,&status);
	else col=1;
	if (fits_read_col(infits,TDOUBLE,col,1,1,cspec->np,&nulval,cspec->wl,
			  &anynul,&status))
	  errormsg("UVES_r1Dspec(): Cannot read wavelength array in file\n\
\t%s",cspec->file);
	/* Read in flux data */
	if (ESO_PHASE3) {
	  if (fits_get_colname(infits,CASEINSEN,"FLUX",card,&col,&status)) {
	    status=0;
	    fits_get_colname(infits,CASEINSEN,"FLUX_EL",card,&col,&status);
	  }
	}
	else col=2;
	if (fits_read_col(infits,TDOUBLE,col,1,1,cspec->np,&nulval,cspec->fl,
			  &anynul,&status))
	  errormsg("UVES_r1Dspec(): Cannot read flux array in file\n\
\t%s",cspec->file);
	/* Read in error data */
	if (ESO_PHASE3) {
	  if (fits_get_colname(infits,CASEINSEN,"ERR",card,&col,&status)) {
	    status=0;
	    fits_get_colname(infits,CASEINSEN,"ERR_EL",card,&col,&status);
	  }
	}
	else col=3;
	if (fits_read_col(infits,TDOUBLE,col,1,1,cspec->np,&nulval,cspec->er,&anynul,&status))
	  errormsg("UVES_r1Dspec(): Cannot read error array in file\n\
\t%s",cspec->file);
	if (naxes>=4 && (!ESO_PHASE3 || UVES_SQUAD_PHASE3)) {
	  /* Read in continuum data */
	  if (UVES_SQUAD_PHASE3) {
	    if (fits_get_colname(infits,CASEINSEN,"CONTINUUM",card,&col,&status))
	      errormsg("UVES_r1Dspec(): Cannot find table column named %s\n\
\tin binary table in FITS file\n\t%s","CONTINUUM",cspec->file);
	  } else col=4;
	  if (fits_read_col(infits,TDOUBLE,col,1,1,cspec->np,&nulval,cspec->co,
			    &anynul,&status))
	    errormsg("UVES_r1Dspec(): Cannot read continuum array in file\n\
\t%s",cspec->file);
	  if (par->nocont==CNTNON) par->nocont=CNTSRD;
	  for (i=0; i<cspec->np; i++) {
	    cspec->no[i]=cspec->fl[i]; cspec->fl[i]*=cspec->co[i];
	    cspec->ne[i]=cspec->er[i]; cspec->er[i]*=cspec->co[i];
	  }
	} else {
	  for (i=0; i<cspec->np; i++) {
	    cspec->no[i]=cspec->fl[i]; cspec->ne[i]=cspec->er[i]; cspec->co[i]=1.0;
	  }
	}
	/* Extract the extra UVES_popler information from UVES_SQUAD Phase 3 files */
	if (UVES_SQUAD_PHASE3) {
	  /* Status */
	  if (fits_get_colname(infits,CASEINSEN,"STATUS",card,&col,&status))
	    errormsg("UVES_r1Dspec(): Cannot find table column named %s\n\
\tin binary table in FITS file\n\t%s","STATUS",cspec->file);
	  if (fits_read_col(infits,TINT,col,1,1,cspec->np,&nulval,cspec->st,
			    &anynul,&status))
	    errormsg("UVES_r1Dspec(): Cannot read status array in file\n\
\t%s",cspec->file);
	  /* Number of pixels before sigma clip */
	  if (fits_get_colname(infits,CASEINSEN,"NPBCLIP",card,&col,&status))
	    errormsg("UVES_r1Dspec(): Cannot find table column named %s\n\
\tin binary table in FITS file\n\t%s","NPBCLIP",cspec->file);
	  if (fits_read_col(infits,TINT,col,1,1,cspec->np,&nulval,cspec->ncb,
			    &anynul,&status))
	    errormsg("UVES_r1Dspec(): Cannot read array of num. contrib.\n\
\tpixels before sigma-clipping in file\n\t%s",cspec->file);
	  /* Number of pixels after sigma clip */
	  if (fits_get_colname(infits,CASEINSEN,"NPACLIP",card,&col,&status))
	    errormsg("UVES_r1Dspec(): Cannot find table column named %s\n\
\tin binary table in FITS file\n\t%s","NPACLIP",cspec->file);
	  if (fits_read_col(infits,TINT,col,1,1,cspec->np,&nulval,cspec->nccb,
			    &anynul,&status))
	    errormsg("UVES_r1Dspec(): Cannot read array of num. contrib.\n\
\tpixels after sigma-clipping in file\n\t%s",cspec->file);
	  /* Chisq. before sigma clip */
	  if (fits_get_colname(infits,CASEINSEN,"CHBCLIP",card,&col,&status))
	    errormsg("UVES_r1Dspec(): Cannot find table column named %s\n\
\tin binary table in FITS file\n\t%s","CHBCLIP",cspec->file);
	  if (fits_read_col(infits,TDOUBLE,col,1,1,cspec->np,&nulval,cspec->csq,
			    &anynul,&status))
	    errormsg("UVES_r1Dspec(): Cannot read array of num. contrib.\n\
\tpixels before sigma-clipping in file\n\t%s",cspec->file);
	  /* Chisq. after sigma clip */
	  if (fits_get_colname(infits,CASEINSEN,"CHACLIP",card,&col,&status))
	    errormsg("UVES_r1Dspec(): Cannot find table column named %s\n\
\tin binary table in FITS file\n\t%s","CHACLIP",cspec->file);
	  if (fits_read_col(infits,TDOUBLE,col,1,1,cspec->np,&nulval,cspec->ccsq,
			    &anynul,&status))
	    errormsg("UVES_r1Dspec(): Cannot read array of num. contrib.\n\
\tpixels after sigma-clipping in file\n\t%s",cspec->file);
	} else {
	  /* Initialize other arrays and check for bad pixels */
	  for (i=0; i<cspec->np; i++) {
	    if (cspec->er[i]==0.0) {
	      cspec->st[i]=RCLIP; cspec->no[i]=1.0; cspec->ncb[i]=cspec->nccb[i]=0;
	      cspec->ef[i]=cspec->er[i]=cspec->nf[i]=cspec->ne[i]=-INFIN;
	    } else {
	      cspec->st[i]=cspec->ncb[i]=cspec->nccb[i]=1;
	      cspec->ef[i]=cspec->er[i]; cspec->nf[i]=cspec->ne[i];
	    }
	    cspec->csq[i]=cspec->ccsq[i]=0.0;
	  }
	}
	/* Make sure there's a reasonable amount of data */
	if (cspec->np<MINUSE)
	  errormsg("UVES_r1Dspec(): Less than MINUSE=%d data points in file\n\
\t%s.\n\tTry reducing MINUSE or zero-pad data file",MINUSE,cspec->file);
	/** Try to determine the dispersion and linear space of the
	    wavelength scale and set wavelengths of pixels edges **/
	/* Currently just a simple test for linearity: compare delta
	   lambda and delta v of first and last pair of wavelength
	   points */
	/* Currently a crude way of setting the edge-wavelength of pixels:
	   take midpoint between successive pixels */
	cdelt=cspec->wl[1]-cspec->wl[0];
	cspec->dwl=cspec->wl[cspec->np-1]-cspec->wl[cspec->np-2];
	if (fabs((cspec->dwl-cdelt)/cspec->dwl)<WLTOL) {
	  par->linear=1; par->disp=cspec->dwl; cspec->dv=0.0;
	  cspec->flwl=cspec->wl[0]-(hdisp=0.5*cspec->dwl);
	  for (i=0; i<cspec->np; i++) cspec->rwl[i]=cspec->wl[i]+hdisp;
	} else {
	  par->linear=0;
	  cdelt=cspec->wl[1]/cspec->wl[0];
	  cspec->dv=cspec->wl[cspec->np-1]/cspec->wl[cspec->np-2];
	  if (fabs((cspec->dv-cdelt)/cspec->dv)>=WLTOL)
	    warnmsg("UVES_r1Dspec(): Wavelength scale in file\n\
\t%s\n\tappears to be neither linear nor log-linear. Assuming log-linear.\n\
\tUVES_popler is currently not equipped to redisperse previously\n\
\tcombined spectra. The wavelength scale of the FITS spectrum written\n\
\tout by UVES_popler will therefore make little sense. However, the ASCII\n\
\tversion of the output file should make sense so I advise using the\n\
\t-dat option. Abandon all hope, ye who enter here.",cspec->file);
	  par->disp=C_C_K*(cspec->dv-1.0);
	  cspec->dv=log10(cspec->dv); cspec->dwl=0.0;
	  cspec->flwl=cspec->wl[0]*pow(10.0,-(hdisp=0.5*cspec->dv));
	  for (i=0; i<cspec->np; i++) cspec->rwl[i]=cspec->wl[i]*pow(10.0,hdisp);
	}
      }
    }
    /* Close the FITS file */
    fits_close_file(infits,&status);
  } else {
    /* File must be an ASCII file. Format is either wl, fl, er or wl,
       nf, ne, co */
    if ((naxes=sscanf(buffer,"%lf %lf %lf",&ddum,&ddum,&ddum))!=3) {
      if ((naxes=sscanf(buffer,"%lf %lf %lf %lf",&ddum,&ddum,&ddum,&ddum))!=4)
	fclose(data_file);
	errormsg("UVES_r1Dspec(): Number of columns in file\n\t%s\n\
\tis %d, but it should be either 3 or 4 with format wl, fl, er\n\
\tor wl norm. fl, norm. er, co",cspec->file);
    }
    /* Find number of lines in file */
    idx++; while (fgets(buffer,LNGSTRLEN,data_file)!=NULL) idx++;
    if (!feof(data_file)) {
      fclose(data_file);
      errormsg("UVES_r1Dspec(): Problem reading line %d in file\n\t%s",idx,
	       cspec->file);
    }
    else { rewind(data_file); cspec->np=idx-1; }
    /* Initialize the combined spectrum but do not set its wavelength scale */
    if (!UVES_init_cspec(cspec,par,1)) {
      fclose(data_file);
      nferrormsg("UVES_r1Dspec(): Error returned from\n\tUVES_init_cspec()");
      return 0;
    }
    /* Read in data */
    for (i=0; i<cspec->np; i++) {
      fgets(buffer,LNGSTRLEN,data_file);
      if (naxes==3) {
	if (sscanf(buffer,"%lf %lf %lf",&(cspec->wl[i]),&(cspec->fl[i]),
		 &(cspec->er[i]))!=3) {
	  fclose(data_file);
	  errormsg("UVES_r1Dspec(): Incorrect format in line %d of file\n\t%s",
		   i+1,cspec->file);
	}
	cspec->no[i]=cspec->fl[i]; cspec->ne[i]=cspec->er[i]; cspec->co[i]=1.0;
      } else if (naxes==4) {
	if (sscanf(buffer,"%lf %lf %lf %lf",&(cspec->wl[i]),&(cspec->no[i]),
		   &(cspec->ne[i]),&(cspec->co[i]))!=4) {
	  fclose(data_file);
	  errormsg("UVES_r1Dspec(): Incorrect format in line %d of file\n\t%s",
		   i+1,cspec->file);
	}
	if (par->nocont==CNTNON) par->nocont=CNTSRD;
	cspec->fl[i]=cspec->no[i]*cspec->co[i];
	cspec->er[i]=cspec->ne[i]*cspec->co[i];
      }
      /* Initialize other arrays */
      if (cspec->er[i]==0.0) {
	cspec->st[i]=RCLIP; cspec->no[i]=1.0; cspec->ncb[i]=cspec->nccb[i]=0;
	cspec->ef[i]=cspec->er[i]=cspec->nf[i]=cspec->ne[i]=-INFIN;
      } else {
	cspec->st[i]=cspec->ncb[i]=cspec->nccb[i]=1;
	cspec->ef[i]=cspec->er[i]; cspec->nf[i]=cspec->ne[i];
      }
      cspec->csq[i]=cspec->ccsq[i]=0.0;
    }
    /* Close the file */
    fclose(data_file);
    /* Make sure there's a reasonable amount of data */
    if (cspec->np<MINUSE)
      errormsg("UVES_r1Dspec(): Less than MINUSE=%d data points in file\n\
\t%s.\n\tTry reducing MINUSE or zero-pad data file",MINUSE,cspec->file);
    /** Try to determine the dispersion and linear space of the
        wavelength scale and set wavelengths of pixels edges **/
    /* Currently just a simple test for linearity: compare delta
       lambda and delta v of first and last pair of wavelength
       points */
    /* Currently a crude way of setting the edge-wavelength of pixels:
       take midpoint between successive pixels */
    cdelt=cspec->wl[1]-cspec->wl[0];
    cspec->dwl=cspec->wl[cspec->np-1]-cspec->wl[cspec->np-2];
    if (fabs((cspec->dwl-cdelt)/cspec->dwl)<WLTOL) {
      par->linear=1; par->disp=cspec->dwl; cspec->dv=0.0;
      cspec->flwl=cspec->wl[0]-(hdisp=0.5*cspec->dwl);
      for (i=0; i<cspec->np; i++) cspec->rwl[i]=cspec->wl[i]+hdisp;
    } else {
      par->linear=0;
      cdelt=cspec->wl[1]/cspec->wl[0];
      cspec->dv=cspec->wl[cspec->np-1]/cspec->wl[cspec->np-2];
      if (fabs((cspec->dv-cdelt)/cspec->dv)>=WLTOL)
	warnmsg("UVES_r1Dspec(): Wavelength scale in file\n\
\t%s\n\tappears to be neither linear nor log-linear. Assuming log-linear.\n\
\tUVES_popler is currently not equipped to redisperse previously\n\
\tcombined spectra. The wavelength scale of the FITS spectrum written\n\
\tout by UVES_popler will therefore make little sense. However, the ASCII\n\
\tversion of the output file should make sense so I advise using the\n\
\t-dat option. Abandon all hope, ye who enter here.",cspec->file);
      par->disp=C_C_K*(cspec->dv-1.0); cspec->dv=log10(cspec->dv);
      cspec->dwl=0.0; 
      cspec->flwl=cspec->wl[0]*pow(10.0,-(hdisp=0.5*cspec->dv));
      for (i=0; i<cspec->np; i++) cspec->rwl[i]=cspec->wl[i]*pow(10.0,hdisp);
    }
  }

  return 1;

}
