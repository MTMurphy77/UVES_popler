/****************************************************************************
* Write out a FITS file containing the raw data (wavelength, flux,
* error etc.) for all echelle orders in a single 1-D spectrum.
****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <fitsio.h>
#include <longnam.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

#define FITS_WKEYD(KEYWORD,COMMENT) \
  if (fits_write_key(outfits,TDOUBLE,KEYWORD,&vald,COMMENT,&status)) \
    errormsg("UVES_wrawFITS(): Cannot write header card %s to file\n\t%s", \
	     KEYWORD,fitsname);
#define FITS_WKEYI(KEYWORD,COMMENT) \
  if (fits_write_key(outfits,TINT,KEYWORD,&vali,COMMENT,&status)) \
    errormsg("UVES_wrawFITS(): Cannot write header card %s to file\n\t%s", \
	     KEYWORD,fitsname);
#define FITS_WKEYS(KEYWORD,COMMENT) \
  if (fits_write_key(outfits,TSTRING,KEYWORD,vals,COMMENT,&status)) \
    errormsg("UVES_wrawFITS(): Cannot write header card %s to file\n\t%s", \
	     KEYWORD,fitsname);

int UVES_wrawFITS(spectrum *spec, params *par) {

  double   crval1=1.0,cd1_1=1.0,cd2_2=1.0,vald=0.0;
  double   *dat=NULL;
  double   version=0.0;
  long     naxes[2]={0,4};
  long     fpixel[2]={1,0};
  int      np=0,vali=0;
  int      crpix1=1,dc_flag=0;
  int      status=0,naxis=2,bitpix=DOUBLE_IMG;
  int      i=0,j=0;
  char     ctype1[FLEN_KEYWORD]="LINEAR\0";
  char     fitsname[LNGSTRLEN]="\0";
  char     clobname[LNGSTRLEN+1]="\0";
  char     comment[FLEN_COMMENT]="\0";
  char     vals[FLEN_VALUE]="\0";
  char     extname[FLEN_COMMENT]="MainFluxFileHeader\0";
  char     *cptr;
  fitsfile *outfits=NULL;
  char     **type=NULL,**form=NULL,**unit=NULL;

  /* Determine output file name based on original file name */
  sprintf(clobname,"!%s",spec->abfile);
  if ((cptr=strrchr(clobname,'.'))==NULL)
    sprintf(clobname,"!%s_up_%03d.fits",spec->abfile,spec->id+1);
  else sprintf(cptr,"_up_%03d.fits",spec->id+1);
  cptr=strchr(clobname,'!')+1; sprintf(fitsname,"%s",cptr);

  /* Determine the maximum number of pixels in any order */
  for (i=0; i<spec->nor; i++) np=(MAX(np,spec->or[i].np));
  /* Allocate memory for temporary data array */
  if ((dat=darray(np))==NULL)
    errormsg("UVES_wrawFITS(): Cannot allocate memory for dat\n\
\tarray of length %d for FITS file\n\t%s",np,fitsname);
  /* Define image dimensions */
  naxes[0]=np; // Other dimensions initialized in variable list

  /* Open specified file */
  if (!access(fitsname,R_OK))
    warnmsg("UVES_wrawFITS(): Overwriting existing file\n\t%s",fitsname);
  else fprintf(stderr,"UVES_wrawFITS(): Opening UVES_popler FITS file\n\
\t%s\n",fitsname);
  if (fits_create_file(&outfits,clobname,&status))
    errormsg("UVES_wrawFITS(): Cannot create FITS file\n\t%s",fitsname);

  /* Create image arrays */
  for (i=0; i<spec->nor; i++) {
    if (fits_create_img(outfits,bitpix,naxis,naxes,&status)) {
errormsg("UVES_wrawFITS(): Cannot create image array %d in\n\
\tFITS file\n\t%s",i+1,fitsname);
    }
  }
  /* Move back to first HDU (not sure why this is needed
     ... "undocumented feature"?) */
  if (fits_movabs_hdu(outfits,1,NULL,&status))
    errormsg("UVES_wrawFITS(): Cannot move to HDU number %d in file\n\
\t%s",1,fitsname);

  /** Write header **/
  /* Leave mark of UVES_popler */
  sprintf(comment,
	  "Output from UVES POst Pipeline Echelle Reduction, by M. T. Murphy,");
  if (fits_write_history(outfits,comment,&status))
    errormsg("UVES_wrawFITS(): Cannot write history keyword\n\tto file %s",
	     fitsname);
  /* Decide which version number to write to output FITS file, either
     the version number in the input UPL file or the version number of
     the running version of UVES_popler */
  version=(par->backvers) ? par->version : VERSION;
  sprintf(comment,"UVES_popler: Version %5.2lf",version);
  if (fits_write_history(outfits,comment,&status))
    errormsg("UVES_wrawFITS(): Cannot write history keyword\n\tto file %s",
	     fitsname);
  sprintf(comment,"%s",WWW);
  if (fits_write_history(outfits,comment,&status))
    errormsg("UVES_wrawFITS(): Cannot write history keyword\n\tto file %s",
	     fitsname);
  /* Date of writing */
  if (fits_write_date(outfits,&status))
    errormsg("UVES_wrawFITS(): Cannot write date to file %s",fitsname);
  /* Number of HDUs(=Number of echelle orders + 1 empty binary table
     with main flux file header) */
  vali=spec->nor+1; FITS_WKEYI("NHDU","Number of HDUs in this file");
  /* Starting pixel */
  vali=crpix1; FITS_WKEYI("CRPIX1","Starting pixel (1-indexed)");
  vald=crval1; FITS_WKEYD("CRVAL1","Starting pixel's dispersion value");
  /* Log-linear flag */
  sprintf(vals,"%s",ctype1);
  FITS_WKEYS("CTYPE1","Log-lin flag: Lin='LINEAR', Log-lin='LOGLIN'");
  vali=dc_flag;
  FITS_WKEYI("DC-FLAG","Log-linear flag: Linear=0, Log-linear=1");
  /* Dispersion */
  vald=cd1_1; FITS_WKEYD("CD1_1","Dispersion of first array");
  vald=cd1_1; FITS_WKEYD("CDELT1","Dispersion of first array");
  vald=cd2_2; FITS_WKEYD("CD2_2","Dispersion of second array");
  vald=cd2_2; FITS_WKEYD("CDELT2","Dispersion of second array");
  /* Object name */
  sprintf(vals,"%s",spec->obj);
  FITS_WKEYS("OBJECT","Astronomical object name");
  /* File name */
  sprintf(vals,"%s",spec->abfile);
  if (strlen(vals)>=FLEN_VALUE)
    errormsg("UVES_wrawFITS(): Length of abbreviated file name,\n\t%s,\n\
\ttoo long (=%d) to be written to output fits file\n\t%s.\n\
\tConsider decreasing file name length and/or increasing length of\n\
\tinternal variable beyond FLEN_VALUE (=%d)",vals,strlen(vals),fitsname,
	     FLEN_VALUE);
  FITS_WKEYS("REDFILE","Reduced flux file name");
  /* Archival file name */
  if (spec->ftype==FTUVES) {
    sprintf(vals,"%s",spec->arfile);
    FITS_WKEYS("ARFILE","Original, archival file name");
  }
  /* File type */
  vali=spec->ftype; FITS_WKEYI("FILETYPE","File type - see web for code");
  /* Number of echelle orders */
  vali=spec->nor; FITS_WKEYI("NECHORD","Number of echelle orders");
  /* Observing date */
  sprintf(vals,"%04d-%02d-%02d",spec->year,spec->month,spec->day);
  FITS_WKEYS("OBSDATE","Observation date (exposure start)");
  /* UT */
  vald=spec->ut; FITS_WKEYD("UT","Universal time (exposure start)");
  /* Julian date */
  vald=spec->jd; FITS_WKEYD("JD","Julian date (exposure start)");
  /* Epoch */
  vald=spec->epoch; FITS_WKEYD("EPOCH","Epoch");
  /* Exposure time */
  vald=spec->etime; FITS_WKEYD("EXPTIME","Exposure time [s]");
  /* Nominal central wavelength */
  vald=spec->cwl; FITS_WKEYD("CENWAV","Nominal central wavelength [nm]");
  /* RA */
  vald=spec->ra; FITS_WKEYD("RA","Right ascension [hrs]");
  /* DEC */
  vald=spec->dec; FITS_WKEYD("DEC","Declination [deg]");
  /* Equinox */
  vald=spec->equ; FITS_WKEYD("EQUINOX","Equinox");
  /* Latitude */
  vald=spec->lat; FITS_WKEYD("OBSLAT","Latitude of observatory [deg]");
  /* Longitude */
  vald=spec->lon; FITS_WKEYD("OBSLON","Longitude of observatory [deg West]");
  /* Altitude */
  vald=spec->alt; FITS_WKEYD("OBSALT","Altitude of observatory [m]");
  /* Heliocentric velocity */
  vald=spec->vhel; FITS_WKEYD("HELVEL","Heliocentric velocity [km/s]");
  /* Velocity shift */
  vald=spec->vshift; FITS_WKEYD("VSHIFT","User-supplied velocity shift [km/s]");
  /* Velocity slope */
  vald=spec->vslope; FITS_WKEYD("VSLOPE","User-supplied velocity slope [km/s/1000Ang]");
  /* Reference wavelength */
  vald=spec->refwav; FITS_WKEYD("REFWAV","User-supplied reference wavelength [Ang]");
  /* Temperature near echelle grating */
  vald=spec->temp; FITS_WKEYD("TEMPGRAT","Temperature near echelle grating [C]");
  /* Atmospheric pressure */
  vald=spec->pres; FITS_WKEYD("ATMOPRES","Atmospheric pressure [hPa]");
  /* Slit width */
  vald=spec->sw; FITS_WKEYD("SLITWID","Slit width [arcsec]");
  /* Reduction pipeline version */
  vali=spec->fvers; FITS_WKEYI("PIPEVERS","UVES Pipeline version (0=MIDAS; 1=CPL)");
  /* UVES_popler-specific descriptions of each array */
  sprintf(vals,"%s","Wavelengths of pixel centres");
  FITS_WKEYS("UP_ARR01","Angstroems");
  sprintf(vals,"%s","Flux spectrum");
  FITS_WKEYS("UP_ARR02","Arbitrary units");
  sprintf(vals,"%s","Error spectrum");
  FITS_WKEYS("UP_ARR03","Arbitrary units; No info. in pixel if <0.0");
  sprintf(vals,"%s","Status");
  FITS_WKEYS("UP_ARR04","1=valid; 0=no data; -ve integer=clipped - see web for code");
  sprintf(vals,"%s","Wavelengths of pixel centres");
  FITS_WKEYS("ARRAY1","Angstroems");
  sprintf(vals,"%s","Flux spectrum");
  FITS_WKEYS("ARRAY2","Arbitrary units");
  sprintf(vals,"%s","Error spectrum");
  FITS_WKEYS("ARRAY3","Arbitrary units; No info. in pixel if <0.0");
  sprintf(vals,"%s","Status");
  FITS_WKEYS("ARRAY4","1=valid; 0=no data; -ve integer=clipped - see web for code");

  /** Fill primary array with first order's data, then fill subsequent
      HDUs with other orders' data **/
  for (i=0; i<spec->nor; i++) {
    /* Move to the correct HDU */
    if (fits_movabs_hdu(outfits,i+1,NULL,&status))
      errormsg("UVES_wrawFITS(): Cannot move to HDU number %d in file\n\
\t%s",i+1,fitsname);
    /* Wavelength */
    fpixel[1]=1;
    if (fits_write_pix(outfits,TDOUBLE,fpixel,spec->or[i].np,spec->or[i].vhwl,
		       &status))
      errormsg("UVES_wrawFITS(): Cannot write wavelength array of length %d\n\
\tfrom order %d of spectrum in file\n\t%s\n\tto FITS file %s",np,i+1,
	       spec->file,fitsname);
    /* Flux */
    fpixel[1]=2;
    if (fits_write_pix(outfits,TDOUBLE,fpixel,spec->or[i].np,spec->or[i].fl,
		       &status))
      errormsg("UVES_wrawFITS(): Cannot write flux array of length %d\n\
\tfrom order %d of spectrum in file\n\t%s\n\tto FITS file %s",np,i+1,
	       spec->file,fitsname);
    /* Error */
    fpixel[1]=3;
    if (fits_write_pix(outfits,TDOUBLE,fpixel,spec->or[i].np,spec->or[i].er,
		       &status))
      errormsg("UVES_wrawFITS(): Cannot write error array of length %d\n\
\tfrom order %d of spectrum in file\n\t%s\n\tto FITS file %s",np,i+1,
	       spec->file,fitsname);
    /* Status */
    fpixel[1]=4;
    for (j=0; j<spec->or[i].np; j++) dat[j]=(double)spec->or[i].st[j];
    for (j=spec->or[i].np; j<np; j++) dat[j]=0.0;
    if (fits_write_pix(outfits,TDOUBLE,fpixel,np,dat,&status))
      errormsg("UVES_wrawFITS(): Cannot write status array of length %d\n\
\tfrom order %d of spectrum in file\n\t%s\n\tto FITS file %s",np,i+1,
	       spec->file,fitsname);
  }

  /* Add the header information from the original flux file to an empty binary table HDU */
  /* Create table */
  if (fits_create_tbl(outfits,BINARY_TBL,0,0,type,form,unit,extname,&status))
    errormsg("UVES_wrawFITS(): Cannot create binary table extension\n\
\tin FITS file\n\t%s",fitsname);
  sprintf(vals,"%s","Header separator 1");
  FITS_WKEYS("HEADSEP1","Main header of original flux file starts below");
  for (i=0; i<spec->nhead_ori; i++) {
    if (fits_write_record(outfits,spec->head_ori[i],&status))
	errormsg("UVES_wrawFITS(): Cannot write header key %d from main\n\
\tflux file %s to raw FITS file\n\t%s",i+1,spec->file,fitsname);
  }

  /* Close file */
  fits_close_file(outfits,&status);

  /* Clean up */
  free(dat);

  return 1;

}
