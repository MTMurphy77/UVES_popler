/****************************************************************************
* Read in a FITS echelle spectrum that has been reduced with PypeIt
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
#include "stdbool.h"

int UVES_r2Dspec_pypeit(spectrum *spec, params *par) {

  double   hour=0.0,min=0.0,sec=0.0;
  double   nulval=0.0,tmpvhel=0.0;
//  double   nrmt[2],onorm=0.0;
  double   dwl=0.0;
  double   *blzfit=NULL;
//  double   *coeff=NULL,*coeffo=NULL;
  bool     badorder=false,longslit=false;
//  long     nrows=0,No=0;
//  long     naxes[9]={0,0,0,0,0,0,0,0,0};
//  int      col=0,naxis=0;
  int      hdutype=0,hdunum=0,status=0,anynul=0;
  int      i=0,j=0,k=0;
  int      thisorder=0,norm=0;
  int      wlcol=0,flcol=0,ercol=0,blzcol=0;
//  int      *ord=NULL;
  char*    slitstr="\0";
  char     pypeline[FLEN_KEYWORD]="Echelle";
  char     comment[FLEN_COMMENT]="\0";
  char     date[FLEN_KEYWORD]="\0",time[FLEN_KEYWORD]="\0";
  char     *cptr;
  fitsfile *infits;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,UVES_replace_envinstr(spec->file),READONLY,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  if (fits_get_hdu_type(infits,&hdutype,&status))
    errormsg("UVES_r2Dspec_pypeit(): Cannot get HDU type for file\n\t%s",spec->file);
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_pypeit(): File not a FITS image: %s",spec->file);

  /* Get the number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec_pypeit(): Cannot find number of HDUs in file\n\
\t%s",spec->file);

  /* Get object name */
  if (par->thar<=1) {
    if (fits_read_key(infits,TSTRING,"TARGET",spec->obj,comment,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","OBJECT",spec->file);
  } else sprintf(spec->obj,"thar_wav");
  /* Alter object name to remove some special characters */
  while ((cptr=strchr(spec->obj,' '))!=NULL) *cptr='_';
  while ((cptr=strchr(spec->obj,'+'))!=NULL) *cptr='p';
  while ((cptr=strchr(spec->obj,'-'))!=NULL) *cptr='m';

  /* Read in entire main header as array of strings */
  if (fits_get_hdrspace(infits,&(spec->nhead_ori),NULL,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read number of keys in main\n\
\theader in FITS file\n\t%s.",spec->file);
  if ((spec->head_ori=cmatrix(spec->nhead_ori,FLEN_CARD))==NULL)
    errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for head_ori\n\
\tmatrix of size %dx%d for file %s\n",spec->nhead_ori,FLEN_CARD,spec->file);
  for (i=0; i<spec->nhead_ori; i++) {
    if (fits_read_record(infits,i+1,spec->head_ori[i],&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read header key %d from\n\
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
    /* Determine the pipeline that was used to reduce the data */
    if (fits_read_key(infits,TSTRING,"PYPELINE",&pypeline,comment,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PYPELINE",spec->file);
    if (strcmp(pypeline,"MultiSlit")==0)
      longslit=true;
    slitstr = longslit ? "spectrum" : "echelle order";

    /* Get modified julian day */
    if (fits_read_key(infits,TDOUBLE,"MJD",&(spec->jd),comment,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MJD",spec->file);
    /* Convert to Julian day */
    spec->jd+=2400000.5;
    /* Get date of observation and convert to year, month and day */
    if (fits_read_key(infits,TSTRING,"DATETIME",date,comment,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DATETIME",spec->file);
    if (strlen(date)>10) {
      /* This means that the UT is appended to the date */
      if (sscanf(date,"%d-%d-%dT%lf:%lf:%lf",&(spec->year),&(spec->month),
		 &(spec->day),&hour,&min,&sec)!=6)
	errormsg("UVES_r2Dspec_pypeit(): Cannot read format of keyword\n\
%s=%s in FITS file\n\t%s","DATETIME",date,spec->file);
      /* Convert hours, mins and secs into the UT */
      spec->ut=sec/3600.0+min/60.0+hour;
    } else {
      /* Only the date should be given */
      if (sscanf(date,"%d-%d-%d",&(spec->year),&(spec->month),&(spec->day))!=3)
	errormsg("UVES_r2Dspec_pypeit(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","DATETIME",date,spec->file);
      /* Convert to hours */
      if (!UVES_hex2dec(time,&(spec->ut)))
	errormsg("UVES_r2Dspec_pypeit(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","UTC",time,spec->file);
    }
    /* Get latitude, longitude and altitude of observatory */
    if (fits_read_key(infits,TDOUBLE,"LAT-OBS",&(spec->lat),comment,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","LAT-OBS",spec->file);
    if (fits_read_key(infits,TDOUBLE,"LON-OBS",&(spec->lon),comment,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","LON-OBS",spec->file);
    if (fits_read_key(infits,TDOUBLE,"ALT-OBS",&(spec->alt),comment,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","ALT-OBS",spec->file);
    /* Get object RA, DEC (both in decimal degrees) and equinox */
    if (fits_read_key(infits,TSTRING,"RA",time,comment,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","RA",spec->file);
    /* Convert RA into hours */
    spec->ra/=15.0;
    if (fits_read_key(infits,TSTRING,"DEC",time,comment,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DEC",spec->file);
    if (fits_read_key(infits,TDOUBLE,"EQUINOX",&(spec->alt),comment,&status)) {
      warnmsg("UVES_r2Dspec_pypeit(): Cannot read value of header card %s\n\
\tfrom FITS file %s.\nAssuming EQUINOX=2000.0","EQUINOX",spec->file);
      spec->equ=2000.0; /* This is the default equinox */
      status=0;
    }
  }

  /* Get exposure time */
  if (fits_read_key(infits,TDOUBLE,"EXPTIME",&(spec->etime),comment,&status))
    errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","EXPTIME",spec->file);

  /* Get number of echelle orders & allocate memory for echelle order array */
  if (fits_read_key(infits,TINT,"NSPEC",&(spec->nor),comment,&status))
    errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NSPEC",spec->file);
  if (!(spec->or=(echorder *)malloc((size_t)(spec->nor*sizeof(echorder)))))
    errormsg("Could not allocate memory for %s array of size %d",
	     slitstr, spec->nor);

  /* For PypeIt files, we now need to advance to the first HDU, where the
     extracted spectra (and some relevant header cards) are stored */
    if (fits_movabs_hdu(infits,2,&hdutype,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot move to second HDU in FITS file\n\
\t%s",spec->file);

  /* Read heliocentric velocity correction */
  /* Note: PypeIt files are assumed to have already been
     corrected to the vacuum-heliocentric frame. However, if
     atmospheric lines are to be masked out by the user, then this
     needs to be done in the observer's frame, so the masked regions
     need to have the same heliocentric correction applied to them as
     the input spectra have */
  if (fits_read_key(infits,TDOUBLE,"VEL_CORR",&tmpvhel,comment,&status))
    errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","VEL_CORR",spec->file);
  //spec->vhel=(tmpvhel-1.0)*299792.458;
  spec->vhel=0; /* Set the correction to zero, as the frames are already in the heliocentric frame */

  /* Find number of pixels to read in for each order */
  if (fits_read_key(infits,TINT,"NAXIS2",&(spec->or[0].np),comment,&status))
    errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS2",spec->file);
  for (i=1; i<spec->nor; i++) spec->or[i].np=spec->or[0].np;

  /* Allocate memory for data arrays and fill wavelength, flux arrays */
  for (i=0; i<spec->nor; i++) {
    /* Move to the HDU corresponding to this order */
    if (fits_movabs_hdu(infits,i+2,&hdutype,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot move to HDU %d in FITS file\n\
\t%s",i+1,spec->file);
    /* Unlike most other instruments/pipelines, PypeIt files have
       an actual wavelength array to be read in, so allocate memory
       for the raw wavelength array first */
    if ((spec->or[i].wl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for raw wavel.\n\
\tarray of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].vhwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for vac. wavel.\n\
\tarray of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].vhrwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for right-edge\n\
\tvac-wavel. array of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,
	       spec->file);
    if ((spec->or[i].fl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for flux\n\
\tarray of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].er=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for error\n\
\tarray of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].res=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for resolution\n\
\tarray of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,spec->file);
    if (par->thar<=1 && !norm) {
      /* Allocate memory for matrix to hold blaze */
      if ((blzfit=darray(spec->or[i].np))==NULL)
        errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for blaze function\n\
\tarray of size %d for file\n\t%s",spec->or[0].np,spec->file);
    }
    if (par->thar==1) {
      if ((spec->or[i].th=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for ThAr\n\
\tflux array of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,spec->file);
      if ((spec->or[i].ter=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for ThAr\n\
\terror array of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,
		 spec->file);
    }
    if ((spec->or[i].st=iarray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for status\n\
\tarray of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,spec->file);
    if (par->thar==1) {
      if ((spec->or[i].tst=iarray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_pypeit(): Cannot allocate memory for ThAr status\n\
\tarray of %s %d of size %d for file\n\t%s",slitstr,i+1,spec->or[i].np,spec->file);
    }

    /* Initialise status array */
    for (j=0; j<spec->or[i].np; j++) spec->or[i].st[j]=1;
    if (par->thar==1) for (j=0; j<spec->or[i].np; j++) spec->or[i].tst[j]=1;
  }

  /* Get number of spectral pixels */
  /*
  if (fits_get_num_rows(infits, &npix, &status))
    errormsg("UVES_r2Dspec_pypeit(): Couldn't determine number of spectral pixels for\n\
\tfile %s",spec->file);*/

  /* PypeIt stores each order in a separate HDU, and each
     HDU contains all wavelength, flux, error information.
     Therefore, we will only loop over the HDUs once, reading in
     all information for each order at a time. */
  for (i=0; i<spec->nor; i++) {
    badorder=false;
    /* Move to the HDU corresponding to this order */
    if (fits_movabs_hdu(infits,i+2,&hdutype,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot move to HDU %d in FITS file\n\
\t%s",i+1,spec->file);
    if (hdutype!=BINARY_TBL)
      errormsg("UVES_r2Dspec_pypeit(): Extension %d not a binary table\n\
\tin file\n\t%s",i+1,spec->file);
    if (longslit) thisorder=i+1;
    else {
      if (fits_read_key(infits,TINT,"HIERARCH ECH_ORDER",&thisorder,comment,&status))
        errormsg("UVES_r2Dspec_pypeit(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","HIERARCH ECH_ORDER",spec->file);
    }
  /* Get the column number for all of the data arrays */
    if (fits_get_colnum(infits, CASEINSEN, "OPT_WAVE", &wlcol, &status)) badorder=true;
    if (fits_get_colnum(infits, CASEINSEN, "OPT_COUNTS", &flcol, &status)) badorder=true;
    if (fits_get_colnum(infits, CASEINSEN, "OPT_COUNTS_SIG", &ercol, &status)) badorder=true;
    if (!norm && fits_get_colnum(infits, CASEINSEN, "OPT_FLAT", &blzcol, &status)) badorder=true;
    /* Check if this is a bad order */
    if (badorder) {
      warnmsg("UVES_r2Dspec_pypeit(): Cannot find optimal extraction of %s %d\n\
\tin FITS file\n\t%s",slitstr,thisorder,spec->file);
        /* Mask all pixels in this order - it is bad */
        spec->or[i].nuse=0;
        for (j=0; j<spec->or[i].np; j++) spec->or[i].st[j]=RCLIP;
        status=0;
        continue;
    }
    /* If we make it to here, then we have found the columns for this order */
    /* Read in flux information */
    if (fits_read_col(infits,TDOUBLE,flcol,1,1,spec->or[i].np,&nulval,
              spec->or[i].fl,&anynul,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read flux array for %s %d\n\
\tin file\n\t%s",slitstr,thisorder,spec->file);
    /* Read in error information */
    if (fits_read_col(infits,TDOUBLE,ercol,1,1,spec->or[i].np,&nulval,
              spec->or[i].er,&anynul,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read error array for %s %d\n\
\tin file\n\t%s",slitstr,thisorder,spec->file);
    /* Any zero values will be assigned -INFIN and marked as clipped during reading */
    for (j=0; j<spec->or[i].np; j++) {
      if (spec->or[i].er[j]<=0.0) {
        spec->or[i].st[j]=RCLIP;
        spec->or[i].er[j]=-INFIN;
      }
    }
    /* Read in vacuum wavelength information */
    if (fits_read_col(infits,TDOUBLE,wlcol,1,1,spec->or[i].np,&nulval,
              spec->or[i].wl,&anynul,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read wavelength array for %s %d\n\
\tin file\n\t%s",slitstr,thisorder,spec->file);
    /* Read in the blaze information */
    if (fits_read_col(infits,TDOUBLE,blzcol,1,1,spec->or[i].np,&nulval,
              blzfit,&anynul,&status))
      errormsg("UVES_r2Dspec_pypeit(): Cannot read wavelength array for %s %d\n\
\tin file\n\t%s",slitstr,thisorder,spec->file);
  /** Normalize by the blaze function if necessary **/
    if (par->thar<=1 && !norm) {
      /* Only operate on useful orders */
      for (j=0; j<spec->or[i].np; j++) {
        spec->or[i].fl[j]/=blzfit[j]; spec->or[i].er[j]/=blzfit[j];
      }
    }
  }
  /* Clean up */
  free(blzfit);

  /* Must treat the special case for PypeIt files where wavelengths are
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

  /* Close flux fits file */
  fits_close_file(infits,&status);

  /* If ThAr information is to be read in, do it now */
  if (par->thar==1) {
    errormsg("UVES_r2Dspec_pypeit(): ThAr information not yet implemented");
  }

  /* Set order identity numbers and parameters */
  for (i=0; i<spec->nor; i++) {
    spec->or[i].id=i; spec->or[i].sid=spec->id;
    /* Find the number of useful pixels in the order */
    if (par->thar<=1) {
      for (spec->or[i].nuse=0,j=0; j<spec->or[i].np; j++) {
	/* Have applied an upper threshold for errors here because
	   PypeIt doesn't apply a negative error bar to regions which
	   are extracted but which are known to contain no flux */
	/* if (spec->or[i].er[j]>DRNDTOL) spec->or[i].nuse++; */
	if (spec->or[i].er[j]>DRNDTOL && spec->or[i].er[j]<5000.0) spec->or[i].nuse++;
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

  /* Read in the polynomial wavelength solutions from the extracted ThAr file.
     The ThAr information in this file is extracted differently to the flux and
     so the wavelength polynomials do not apply to the flux and error arrays
     above. Indeed, above we read in the wavelength scale directly, not via
     polynomials like below */
  if (par->thar==2) {
    errormsg("UVES_r2Dspec_pypeit(): ThAr information not yet implemented");
  }

  /* Set median resolution to zero for lack of any other information */
  spec->arcfwhm=0.0;

  /* Make sure this spectrum will be combined into combined spectrum later */
  spec->comb=1;
  return 1;

}
