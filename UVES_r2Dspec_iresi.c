/****************************************************************************
* Read in a IRAF-reduced Keck/ESI 2D FITS echelle spectrum and error array
****************************************************************************/

#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <longnam.h>
#include "UVES_popler.h"
#include "astron.h"
#include "memory.h"
#include "error.h"

int UVES_r2Dspec_iresi(spectrum *spec, params *par) {

  double   hour=0.0,min=0.0,sec=0.0;
  double   nulval=0.0;
  long     naxes[9]={0,0,0,0,0,0,0,0,0};
  int      npix=0,nspace=0;
  int      hdutype=0,hdunum=0,status=0,bitpix=0,first=1,naxis=0,anynul=0;
  int      i=0,j=0,k=0;
  char     comment[FLEN_COMMENT]="\0";
  char     date[FLEN_KEYWORD]="\0",time[FLEN_KEYWORD]="\0";
  char     card[FLEN_CARD]="\0";
  char     form1[LNGSTRLEN]="\0";
  char     dummy[FLEN_KEYWORD]="\0",dummy2[FLEN_KEYWORD]="\0";
  char     buffer[VVVHUGESTRLEN]="\0";
  char     *cptr;
  fitsfile *infits;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,UVES_replace_envinstr(spec->file),READONLY,&status))
      errormsg("UVES_r2Dspec_iresi(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  if (fits_get_hdu_type(infits,&hdutype,&status))
    errormsg("UVES_r2Dspec_iresi(): Cannot get HDU type for file\n\t%s",spec->file);
  if (hdutype!=IMAGE_HDU)
    errormsg("UVES_r2Dspec_iresi(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("UVES_r2Dspec_iresi(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (hdunum>1)
    errormsg("UVES_r2Dspec_iresi(): Number of HDUs is %d instead\n\
\tof %d (at most) in file\n\t%s",hdunum,1,spec->file);

  /* Get object name */
  if (par->thar<=1) {
    if (fits_read_key(infits,TSTRING,"OBJECT",spec->obj,comment,&status))
      errormsg("UVES_r2Dspec_iresi(): Cannot read value of header card \n\
\t%s from FITS file\n\t%s.","OBJECT",spec->file);
  }
  else sprintf(spec->obj,"thar_wav");
  /* Alter object name to remove some special characters */
  while ((cptr=strchr(spec->obj,' '))!=NULL) *cptr='_';
  while ((cptr=strchr(spec->obj,'+'))!=NULL) *cptr='p';
  while ((cptr=strchr(spec->obj,'-'))!=NULL) *cptr='m';

  /* Read in entire main header as array of strings */
  if (fits_get_hdrspace(infits,&(spec->nhead_ori),NULL,&status))
      errormsg("UVES_r2Dspec_iresi(): Cannot read number of keys in main\n\
\theader in FITS file\n\t%s.",spec->file);
  if ((spec->head_ori=cmatrix(spec->nhead_ori,FLEN_CARD))==NULL)
    errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for head_ori\n\
\tmatrix of size %dx%d for file %s\n",spec->nhead_ori,FLEN_CARD,spec->file);
  for (i=0; i<spec->nhead_ori; i++) {
    if (fits_read_record(infits,i+1,spec->head_ori[i],&status))
      errormsg("UVES_r2Dspec_iresi(): Cannot read header key %d from\n\
\tmain header of file %s\n",i+1,spec->file);
  }

  /* Read in archival filename */
  sprintf(spec->arfile,"%s",spec->abfile);

  /* Determine the slit-width, CCD binning and the temperature */
  spec->temp=spec->sw=0.0;
  if (fits_read_key(infits,TSTRING,"BINNING",card,comment,&status))
    errormsg("UVES_r2Dspec_iresi(): Cannot read value of header card\n\
%s from FITS file\n\\t%s.","BINNING",spec->file);
  if (sscanf(card,"%d, %d",&(spec->biny),&(spec->binx))!=2)
    errormsg("UVES_r2Dspec_iresi(): Do not understand header card %s\n\
\t(='%s') from FITS file\n\\t%s.","BINNING",card,spec->file);
  if (par->thar<=1) {
    /* Get date of observation and convert to year, month and day */
    if (fits_read_key(infits,TSTRING,"DATE-OBS",date,comment,&status))
      errormsg("UVES_r2Dspec_iresi(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DATE-OBS",spec->file);
    if (sscanf(date,"%d-%d-%d",&(spec->year),&(spec->month),&(spec->day))!=3)
      errormsg("UVES_r2Dspec_iresi(): Cannot read format of keyword\n\
\t%s=%s in FITS file\n\t%s","DATE-OBS",date,spec->file);
    /* Get universal time from UT */ 
    if (fits_read_key(infits,TSTRING,"UT",card,comment,&status))
      errormsg("UVES_r2Dspec_iresi(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","UT",spec->file);
    /* Convert to hours */
    if (sscanf(card,"%lf:%lf:%lf",&hour,&min,&sec)!=3)
      errormsg("UVES_r2Dspec_iresi(): Do not understand header card %s\n\
\t(=%s) from FITS file\n\t%s.","UT",card,spec->file);
    spec->ut=hour+min/60.0+sec/3600.0;
    /* Convert date+time to Julian day */
    spec->jd=ast_date2jd(spec->year,spec->month,spec->day,spec->ut);
    /* Get latitude, longitude and altitude of observatory */
    /* At the moment, these are hard-coded assuming that we really are
       dealing with Keck/ESI data */
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
  if (fits_read_key(infits,TDOUBLE,"EXPOSURE",&(spec->etime),comment,&status))
    errormsg("UVES_r2Dspec_iresi(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","EXPOSURE",spec->file);

  /* Get number of echelle orders & allocate memory for echelle order array */
  if (fits_read_key(infits,TINT,"NAXIS2",&(spec->nor),comment,&status))
    errormsg("UVES_r2Dspec_iresi(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS2",spec->file);
  if (!(spec->or=(echorder *)malloc((size_t)(spec->nor*sizeof(echorder)))))
    errormsg("Could not allocate memory for echelle order array of size %d",
	     spec->nor);

  /* Find number of pixels to read in for each order */
  if (fits_read_key(infits,TINT,"NAXIS1",&(spec->or[0].np),comment,&status))
    errormsg("UVES_r2Dspec_iresi(): Cannot read value of header card\n\
\t%s from FITS file\n\t%s.","NAXIS1",spec->file);
  for (i=1; i<spec->nor; i++) spec->or[i].np=spec->or[0].np;

  /* Allocate memory for data arrays and fill wavelength, flux arrays */
  for (i=0; i<spec->nor; i++) {
    if ((spec->or[i].vhwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for vac. wavel.\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].vhrwl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for right-edge\n\
\tvac-wavel. array of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,
	       spec->file);
    if ((spec->or[i].fl=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for flux\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].er=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for error\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if ((spec->or[i].res=darray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for resolution\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if (par->thar==1) {
      if ((spec->or[i].th=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for ThAr flux\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
      if ((spec->or[i].ter=darray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for ThAr error\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    }
    if ((spec->or[i].st=iarray(spec->or[i].np))==NULL)
      errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for status\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    if (par->thar==1) {
      if ((spec->or[i].tst=iarray(spec->or[i].np))==NULL)
	errormsg("UVES_r2Dspec_iresi(): Cannot allocate memory for ThAr status\n\
\tarray of order %d of size %d for file\n\t%s",i+1,spec->or[i].np,spec->file);
    }
    /* Initialise status array */
    for (j=0; j<spec->or[i].np; j++) spec->or[i].st[j]=1;
    if (par->thar==1) for (j=0; j<spec->or[i].np; j++) spec->or[i].tst[j]=1;
  }

  /* Get image dimensions */
  if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
    errormsg("UVES_r2Dspec_iresi(): Couldn't get image dimensions for\n\
\tfile %s",spec->file);
  if (!naxis) errormsg("UVES_r2Dspec_iresi(): Couldn't find input image extension\n\
\tfor file %s",spec->file);
  npix=naxes[0];
  /* Read in flux and error information */
  for (i=0,first=1; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].fl,&anynul,&status))
      errormsg("UVES_r2Dspec_iresi(): Cannot read flux array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }
  first=3*spec->nor*npix;
  for (i=0; i<spec->nor; i++,first+=npix) {
    if (fits_read_img(infits,TDOUBLE,first,spec->or[i].np,&nulval,
		      spec->or[i].er,&anynul,&status))
      errormsg("UVES_r2Dspec_iresi(): Cannot read error array for order\n\
\t%d in file\n\t%s",i+1,spec->file);
  }
  /* Read in wavelength solution information */
  /* At the moment it is assumed that each order has been given its
     own linear wavelength dispersion and that this is stored in the
     WAT2_* header cards in the specific format below. Only the
     starting wavelength and linear dispersion are read for each
     order */
  i=1; sprintf(dummy,"WAT2_%03d",i); strcpy(buffer,"\0");
  while (!fits_read_key(infits,TSTRING,dummy,dummy2,comment,&status)) {
    if (strlen(buffer)+strlen(dummy2)>=VVVHUGESTRLEN)
      errormsg("UVES_r2Dspec_iresi(): Cummulative length of WAT2_*\n\
\theader keywords too long in file\n\t%s",spec->file);
    nspace=68-strlen(dummy2);
    strcat(buffer,dummy2); for (j=0; j<nspace; j++) strcat(buffer," ");
    i++; sprintf(dummy,"WAT2_%03d",i);
  }
  /* Break down long string into wavelength coefficients for each order */
  sprintf(form1,"%s","%*d %*d %*d %lf %lf %*d %*lf %*lf %*lf");
  for (i=0; i<spec->nor; i++) {
    spec->or[i].nwpol=2;
    sprintf(dummy,"spec%d = \"",i+1);
    if ((cptr=strstr(buffer,dummy))==NULL)
      errormsg("UVES_r2Dspec_iresi(): Cannot find wavelength coefficients\n\
\tfor order %d in file\n\t%s",i+1,spec->file);
    cptr+=strlen(dummy);
    if (sscanf(cptr,form1,&(spec->or[i].wpol[0]),&(spec->or[i].wpol[1]))!=2)
      errormsg("UVES_r2Dspec_iresi(): Cannot read wavelength\n\
\tcoefficients for order %d in file\n\t%s",i+1,spec->file);
    /* fprintf(stdout,"%4d %lf %lf\n",i+1,spec->or[i].wpol[0],spec->or[i].wpol[1]); */
    for (j=spec->or[i].nwpol; j<NWPOL; j++) spec->or[i].wpol[j]=0.0;
  }

  /* Close flux fits file */
  fits_close_file(infits,&status);

  /* If ThAr information is to be read in, do it now */
  if (par->thar==1)
    errormsg("UVES_r2Dspec_iresi(): Do not yet know how to read in ThAr\n\
\tinformation for IRAF-ESI files");

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
      free(spec->or[i].vhwl); free(spec->or[i].vhrwl); free(spec->or[i].fl);
      free(spec->or[i].er); free(spec->or[i].res); free(spec->or[i].st);
      if (par->thar==1) {
	free(spec->or[i].th); free(spec->or[i].ter); free(spec->or[i].tst);
      }
    }
  }

  /* Make sure this spectrum will be combined into combined spectrum later */
  spec->comb=1;

  return 1;

}
