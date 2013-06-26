/****************************************************************************
* Write out a UVES_popler FITS file containing the 1-D spectrum and,
* if required, the redispersed ThAr spectra.
****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fitsio.h>
#include <longnam.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

int UVES_wFITSfile(spectrum *spec, int nspec, cspectrum *cspec, params *par) {

  double   crval1=0.0,cd1_1=0.0,cd2_2=1.0;
  double   disp=0.0;
  double   *dat=NULL;
  double   dtbl[nspec];
  double   version=0.0;
  long     naxes[2]={0,0};
  long     fpixel[2]={0,0};
  int      crpix1=1,dc_flag=0;
  int      status=0,naxis=2,bitpix=DOUBLE_IMG;
  int      itbl[nspec];
  int      i=0,j=0;
  char     ctype1[FLEN_KEYWORD]="\0";
  char     fitsname[LNGSTRLEN]="\0";
  char     comment[FLEN_COMMENT]="\0";
  char     extname[FLEN_COMMENT]="ExposureInformation\0";
  char *type[TABNCOL]={"FileName","ArFiName","FileType","ObjName", \
		       "DATE-OBS","JD","Epoch","UT","ExpTime","CenWav","RA",\
		       "Dec","Latitude","Longitud","Altitude","HelioVel","Temp",\
		       "Pressure","ThFiName","ThArFiNm","ThArJD",\
		       "MEDFWHM","ThArTemp","ThArPress"};
  char *form[TABNCOL]={"1A","1A","1I","1A","1A","1D","1D","1D","1D","1D",\
		       "1D","1D","1D","1D","1D","1D","1D","1D","1A","1A",\
		       "1D","1D","1D","1D"};
  char *unit[TABNCOL]={" "," "," "," ",\
		       "UT date at exposure start [YYY-MM-DD]",\
		       "Julian date at exp. start [days]",\
		       "Epoch at exposure start [yrs]",\
		       "UT at exposure start [hrs]","Exposure time [s]",\
		       "Nominal central wavelength [nm]","RA of telescope [hrs]",\
		       "Declination of telescope [deg]",\
		       "Latitude of observatory [deg]",\
		       "Longitude: positive=West [deg]",\
		       "Altitude of observatory [m]",\
		       "Heliocentric velocity [km/s]",\
		       "Temperature at grating [deg C or K]",\
		       "Atmospheric pressure [hPa]"," "," ",\
		       "Julian date for ThAr [days]",\
		       "Median FWHM of ThAr lines [km/s]",\
		       "Temp. at grat. for ThAr [deg C or K]",\
		       "Atmospheric pressure for ThAr [hPa]"};
  char     **stbl;
  char     *cptr;
  fitsfile *outfits=NULL;

  /* Allocate memory for temporary data array */
  if ((dat=darray(cspec->np))==NULL)
    errormsg("UVES_wFITSfile(): Cannot allocate memory for dat\n\
\tarray of length %d",cspec->np);

  /* Open specified file */
  if (!access(cspec->FITSfile,R_OK))
    warnmsg("Overwriting existing file %s",cspec->FITSfile);
  else fprintf(stderr,"Opening UVES_popler FITS file %s\n",cspec->FITSfile);
  sprintf(fitsname,"!%s",cspec->FITSfile);
  if (fits_create_file(&outfits,fitsname,&status))
    errormsg("UVES_wFITSfile(): Cannot create FITS file %s",cspec->FITSfile);
  /* Define image dimensions */
  naxes[0]=cspec->np; naxes[1]=9;
  /* Create primary array */
  if (fits_create_img(outfits,bitpix,naxis,naxes,&status))
    errormsg("UVES_wFITSfile(): Cannot create primary array in\n\
\tFITS file\n\t%s",cspec->FITSfile);
  /** Fill primary array **/
  /* Normalized Flux */
  fpixel[0]=1; fpixel[1]=1;
  if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],cspec->no,&status))
    errormsg("UVES_wFITSfile(): Cannot write flux array of length %d\n\
\tto FITS file %s",naxes[0],cspec->FITSfile);
  /* Normalized Error */
  fpixel[0]=1; fpixel[1]=2;
  /***** TEMP!!!!!!! *****/
  /* for (i=0; i<cspec->np; i++) dat[i]=(cspec->ne[i]>0.0) ? cspec->ne[i] : 1000.0; */
  /* if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],dat,&status)) */
  if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],cspec->ne,&status))
    errormsg("UVES_wFITSfile(): Cannot write error array of length %d\n\
\tto FITS file %s",naxes[0],cspec->FITSfile);
  /* Normalized Fluctuation */
  fpixel[0]=1; fpixel[1]=3;
  if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],cspec->nf,&status))
    errormsg("UVES_wFITSfile(): Cannot write error array of length %d\n\
\tto FITS file %s",naxes[0],cspec->FITSfile);
  /* Continuum */
  fpixel[0]=1; fpixel[1]=4;
  if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],cspec->co,&status))
    errormsg("UVES_wFITSfile(): Cannot write cont. array of length %d\n\
\tto FITS file %s",naxes[0],cspec->FITSfile);
  /* Status */
  fpixel[0]=1; fpixel[1]=5;
  for (i=0; i<cspec->np; i++) dat[i]=(double)cspec->st[i];
  if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],dat,&status))
    errormsg("UVES_wFITSfile(): Cannot write status array of length %d\n\
\tto FITS file %s",naxes[0],cspec->FITSfile);
  /* Number of pixels before clip */
  fpixel[0]=1; fpixel[1]=6;
  for (i=0; i<cspec->np; i++) dat[i]=(double)cspec->ncb[i];
  if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],dat,&status))
    errormsg("UVES_wFITSfile(): Cannot write ncb array of length %d\n\
\tto FITS file %s",naxes[0],cspec->FITSfile);
  /* Number of pixels after clip */
  fpixel[0]=1; fpixel[1]=7;
  for (i=0; i<cspec->np; i++) dat[i]=(double)cspec->nccb[i];
  if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],dat,&status))
    errormsg("UVES_wFITSfile(): Cannot write nccb array of length %d\n\
\tto FITS file %s",naxes[0],cspec->FITSfile);
  /* Chisquared before clip */
  fpixel[0]=1; fpixel[1]=8;
  if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],cspec->csq,&status))
    errormsg("UVES_wFITSfile(): Cannot write csq array of length %d\n\
\tto FITS file %s",naxes[0],cspec->FITSfile);
  /* Chisquared after clip */
  fpixel[0]=1; fpixel[1]=9;
  if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],cspec->ccsq,&status))
    errormsg("UVES_wFITSfile(): Cannot write ccsq array of length %d\n\
\tto FITS file %s",naxes[0],cspec->FITSfile);
  /* Leave mark of UVES_popler */
  sprintf(comment,
	  "Output from UVES POst Pipeline Echelle Reduction, by M. T. Murphy,");
  if (fits_write_history(outfits,comment,&status))
    errormsg("UVES_wFITSfile(): Cannot write history keyword about origin\n\
\ttto file %s",cspec->FITSfile);
  /* Decide which version number to write to output FITS file, either
     the version number in the input UPL file or the version number of
     the running version of UVES_popler */
  version=(par->backvers) ? par->version : VERSION;
  sprintf(comment,"UVES_popler: Version %5.2lf",version);
  if (fits_write_history(outfits,comment,&status))
    errormsg("UVES_wFITSfile(): Cannot write UVES_popler version\n\tto file %s",
	     cspec->FITSfile);
  /*
  sprintf(comment,"\t%s",WWW);
  if (fits_write_history(outfits,comment,&status))
    errormsg("UVES_wFITSfile(): Cannot write UVES_popler website\n\tto file %s",
	     cspec->FITSfile);
  */
  /* Date of writing */
  if (fits_write_date(outfits,&status))
    errormsg("UVES_wFITSfile(): Cannot write date to file %s",cspec->FITSfile);
  /* Object name */
  sprintf(comment,"Object name");
  if (fits_write_key(outfits,TSTRING,"OBJECT",cspec->obj,comment,&status))
    errormsg("UVES_wFITSfile(): Cannot write object name to file\n\t%s",
	     cspec->FITSfile);
  /* Write header keywords for wavelength information */
  if (par->linear) {
    crval1=cspec->wl[0]; cd1_1=cspec->dwl; dc_flag=0; sprintf(ctype1,"LINEAR");
  } else {
    crval1=log10(cspec->wl[0]); cd1_1=cspec->dv; dc_flag=1; sprintf(ctype1,"LOGLIN");
  }
  /* Starting pixel */
  sprintf(comment,"Starting pixel (1-indexed)");
  if (fits_write_key(outfits,TINT,"CRPIX1",&crpix1,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add CRPIX1 keyword to file\n\t%s",
	     cspec->FITSfile);
  /* Log-linear flag */
  sprintf(comment,"Log-lin flag: Lin='LINEAR', Log-lin='LOGLIN'");
  if (fits_write_key(outfits,TSTRING,"CTYPE1",ctype1,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add CTYPE1 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"Log-linear flag: Linear=0, Log-linear=1");
  if (fits_write_key(outfits,TINT,"DC-FLAG",&dc_flag,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add DC_FLAG keyword to file\n\t%s",
	     cspec->FITSfile);
  /* Starting wavelength */
  if (par->linear) sprintf(comment,"Center wavelength of first pixel");
  else sprintf(comment,"Center wavelength (log10) of first pixel");
  if (fits_write_key(outfits,TDOUBLE,"CRVAL1",&crval1,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add CRVAL1 keyword to file\n\t%s",
	     cspec->FITSfile);
  /* Dispersion */
  if (par->linear) sprintf(comment,"Dispersion: w[i]=w[0]+([i+1]-CRPIX1)*CD1_1");
  else sprintf(comment,"Dispersion: w[i]=wl[0]*10^{([i+1]-CRPIX1)*CD1_1}");
  if (fits_write_key(outfits,TDOUBLE,"CD1_1",&cd1_1,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add CD1_1 keyword to file\n\t%s",
	     cspec->FITSfile);
  if (fits_write_key(outfits,TDOUBLE,"CDELT1",&cd1_1,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add CDELT1 keyword to file\n\t%s",
	     cspec->FITSfile);
  if (fits_write_key(outfits,TDOUBLE,"CD2_2",&cd2_2,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add CD2_2 keyword to file\n\t%s",
	     cspec->FITSfile);
  if (fits_write_key(outfits,TDOUBLE,"CDELT2",&cd2_2,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add CDELT2 keyword to file\n\t%s",
	     cspec->FITSfile);
  /* UVES_popler-specific values */
  sprintf(comment,"Center wavelength of first pixel [Angstroms]");
  if (fits_write_key(outfits,TDOUBLE,"UP_WLSRT",&cspec->wl[0],comment,&status))
    errormsg("UVES_wFITSfile(): Could not add UP_WLSRT keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"Center wavelength of last pixel [Angstroms]");
  if (fits_write_key(outfits,TDOUBLE,"UP_WLEND",&cspec->wl[cspec->np-1],comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add UP_WLEND keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"Emission redshift");
  if (fits_write_key(outfits,TDOUBLE,"UP_ZEM",&par->zem,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add UP_ZEM keyword to file\n\t%s",
	     cspec->FITSfile);
  if (par->linear) { disp=cd1_1; sprintf(comment,"Dispersion [Angstroms]"); }
  else { disp=par->disp; sprintf(comment,"Dispersion [km/s]"); }
  if (fits_write_key(outfits,TDOUBLE,"UP_DISP",&disp,comment,&status))
    errormsg("UVES_wFITSfile(): Could not add UP_DISP keyword to file\n\t%s",
	     cspec->FITSfile);
  /* Write descriptive keywords for each array */
  sprintf(comment,"No units");
  if (fits_write_key(outfits,TSTRING,"UP_ARR01","Normalized flux spectrum",
		     comment,&status))
    errormsg("UVES_wFITSfile(): Could not add UP_ARR01 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"No units; No info. in pixel if <0.0");
  if (fits_write_key(outfits,TSTRING,"UP_ARR02","Normalized error",
		     comment,&status))
    errormsg("UVES_wFITSfile(): Could not add UP_ARR02 keyword to file\n\t%s",
	     cspec->FITSfile);
  if (fits_write_key(outfits,TSTRING,"UP_ARR03","Norm. expected fluctuation",
		     comment,&status))
    errormsg("UVES_wFITSfile(): Could not add UP_ARR03 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"Arbitrary units");
  if (fits_write_key(outfits,TSTRING,"UP_ARR04",
		     "Continuum by which data are normalized",comment,&status))
    errormsg("UVES_wFITSfile(): Could not add UP_ARR04 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"1=valid; -ve integer=clipped - see web for code");
  if (fits_write_key(outfits,TSTRING,"UP_ARR05","Status",comment,&status))
    errormsg("UVES_wFITSfile(): Could not add UP_ARR05 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"# pix contrib. to comb. pix. before sigma-clip");
  if (fits_write_key(outfits,TSTRING,"UP_ARR06","NPIX before clip",comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add UP_ARR06 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"# pix contrib. to comb. pix. after sigma-clip");
  if (fits_write_key(outfits,TSTRING,"UP_ARR07","NPIX after clip",comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add UP_ARR07 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"chisq. about weighted mean before sigma-clip");
  if (fits_write_key(outfits,TSTRING,"UP_ARR08","chisq. before clip",comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add UP_ARR08 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"chisq. about weighted mean after sigma-clip");
  if (fits_write_key(outfits,TSTRING,"UP_ARR09","chisq. after clip",comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add UP_ARR09 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"No units");
  if (fits_write_key(outfits,TSTRING,"ARRAY1","Normalized flux spectrum",comment,&status))
    errormsg("UVES_wFITSfile(): Could not add ARRAY1 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"No units; No info. in pixel if <0.0");
  if (fits_write_key(outfits,TSTRING,"ARRAY2","Normalized error",comment,&status))
    errormsg("UVES_wFITSfile(): Could not add ARRAY2 keyword to file\n\t%s",
	     cspec->FITSfile);
  if (fits_write_key(outfits,TSTRING,"ARRAY3","Norm. expected fluctuation",comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add ARRAY3 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"Arbitrary units");
  if (fits_write_key(outfits,TSTRING,"ARRAY4",
		     "Continuum by which data are normalized",comment,&status))
    errormsg("UVES_wFITSfile(): Could not add ARRAY4 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"1=valid; -ve integer=clipped - see web for code");
  if (fits_write_key(outfits,TSTRING,"ARRAY5","Status",comment,&status))
    errormsg("UVES_wFITSfile(): Could not add ARRAY5 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"# pix contrib. to comb. pix. before sigma-clip");
  if (fits_write_key(outfits,TSTRING,"ARRAY6","NPIX before clip",comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add ARRAY6 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"# pix contrib. to comb. pix. after sigma-clip");
  if (fits_write_key(outfits,TSTRING,"ARRAY7","NPIX after clip",comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add ARRAY7 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"chisq. about weighted mean before sigma-clip");
  if (fits_write_key(outfits,TSTRING,"ARRAY8","chisq. before clip",comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add ARRAY8 keyword to file\n\t%s",
	     cspec->FITSfile);
  sprintf(comment,"chisq. about weighted mean after sigma-clip");
  if (fits_write_key(outfits,TSTRING,"ARRAY9","chisq. after clip",comment,
		     &status))
    errormsg("UVES_wFITSfile(): Could not add ARRAY9 keyword to file\n\t%s",
	     cspec->FITSfile);

  /** Create binary table HDU containing parameters from each exposure **/
  /* Only do this if we have contributing spectra, i.e. not for
     already-combined spectra */
  if (nspec) {
    /* Create table */
    if (fits_create_tbl(outfits,BINARY_TBL,0,TABNCOL,type,form,unit,extname,
			&status))
      errormsg("UVES_wFITSfile(): Cannot create binary table extension\n\
\tin FITS file\n\t%s",cspec->FITSfile);
    /* Allocate memory for string matrix */
    if ((stbl=cmatrix(nspec,FLEN_KEYWORD))==NULL)
      errormsg("UVES_wFITSfile(): Cannot allocate memory for stble\n\
\tarray of length %dx%d",nspec,FLEN_KEYWORD);
    /* Fill table rows */
    j=0;
    /* File name */
    for (i=0; i<nspec; i++) {
      if (strlen(spec[i].abfile)<FLEN_KEYWORD)
	sprintf(stbl[i],"%s",spec[i].abfile);
      else errormsg("UVES_wFITSfile(): Length of abbreviated file name,\n\t%s,\n\
\ttoo long (=%d) to be written to output fits file\n\t%s.\n		\
\tConsider decreasing file name length and/or increasing length of\n	\
\tinternal variable beyond FLEN_KEYWORD (=%d)",spec[i].abfile,
		    strlen(spec[i].abfile),cspec->FITSfile,FLEN_KEYWORD);
    }
    if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Archival file name */
    for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].arfile);
    if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* File type */
    for (i=0; i<nspec; i++) itbl[i]=spec[i].ftype;
    if (fits_write_col(outfits,TINT,j+1,1,1,nspec,&itbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Object name */
    for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].obj);
    if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Date */
    for (i=0; i<nspec; i++)
      sprintf(stbl[i],"%04d-%02d-%02d",spec[i].year,spec[i].month,spec[i].day);
    if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Julian date */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].jd;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Epoch */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].epoch;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* UT */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].ut;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Exposure time */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].etime;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Nominal central wavelength */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].cwl;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* RA */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].ra;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* DEC */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].dec;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Latitude */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].lat;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Longitude */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].lon;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Altitude */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].alt;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Heliocentric velocity */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].vhel;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Temperature */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].temp;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Pressure */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].pres;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* ThAr file name */
    for (i=0; i<nspec; i++) {
      if (strlen(spec[i].abthfile)<FLEN_KEYWORD)
	sprintf(stbl[i],"%s",spec[i].abthfile);
      else errormsg("UVES_wFITSfile(): Length of abbrev. ThAr file name,\n\t%s,\n\
\ttoo long (=%d) to be written to output fits file\n\t%s.\n		\
\tConsider decreasing file name length and/or increasing length of\n	\
\tinternal variable beyond FLEN_KEYWORD (=%d)",spec[i].abthfile,
		    strlen(spec[i].abfile),cspec->FITSfile,FLEN_KEYWORD);
    }
    if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status))
    errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* ThAr archival file name */
    for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].tharfile);
    if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* ThAr Julian date */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].thjd;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* Median FWHM */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].arcfwhm;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* ThAr temperature */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].thtemp;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    j++;
    /* ThAr pressure */
    for (i=0; i<nspec; i++) dtbl[i]=spec[i].thpres;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,&dtbl,&status))
      errormsg("UVES_wFITSfile(): Cannot add column '%s' to \n\
\tbinary table extension in FITS file\n\t%s",type[j],cspec->FITSfile);
    /* Clean up */
    free(*stbl); free(stbl);
  }
  /* Close file */
  fits_close_file(outfits,&status);

  /* If required, write out a FITS file containing all the
     redispersed, merged ThAr spectra */
  if (par->thar==1) {
    /* Define name of ThAr output file */
    sprintf(fitsname,"%s",cspec->FITSfile); cptr=strstr(fitsname,".fits");
    sprintf(cptr,"_thar.fits");
    /* Open specified file */
    if (!access(fitsname,R_OK)) warnmsg("Overwriting existing file %s",fitsname);
    else fprintf(stderr,"Opening UVES_popler ThAr FITS file %s\n",fitsname);
    sprintf(fitsname,"!%s",cspec->FITSfile); cptr=strstr(fitsname,".fits");
    sprintf(cptr,"_thar.fits"); cptr=fitsname+1;
    if (fits_create_file(&outfits,fitsname,&status))
      errormsg("UVES_wFITSfile(): Cannot create FITS file %s",cptr);
    for (j=0; j<nspec; j++) {
      /* Define image dimensions */
      naxes[0]=spec[j].th.np; naxes[1]=5;
      /* Create primary array */
      if (fits_create_img(outfits,bitpix,naxis,naxes,&status))
	errormsg("UVES_wFITSfile(): Cannot create primary array in FITS\n\t%s",
		 cptr);
      /* Normalized ThAr Flux */
      fpixel[0]=1; fpixel[1]=1;
      if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],spec[j].th.fl,&status))
	errormsg("UVES_wFITSfile(): Cannot write ThAr flux array of length %d\n\
\tto FITS file %s",naxes[0],cptr);
      /* Normalized ThAr Error */
      fpixel[0]=1; fpixel[1]=2;
      if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],spec[j].th.er,&status))
	errormsg("UVES_wFITSfile(): Cannot write ThAr error array of length %d\n\
\tto FITS file %s",naxes[0],cptr);
      /* ThAr Status */
      fpixel[0]=1; fpixel[1]=3;
      for (i=0; i<spec[j].th.np; i++) dat[i]=(double)spec[j].th.st[i];
      if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],dat,&status))
	errormsg("UVES_wFITSfile(): Cannot write ThAr status array of length %d\n\
\tto FITS file %s",naxes[0],cptr);
      /* Normalized Error */
      fpixel[0]=1; fpixel[1]=4;
      if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],
			 &(cspec->ne[spec[j].th.csidx]),&status))
	errormsg("UVES_wFITSfile(): Cannot write error array of length %d\n\
\tto FITS file %s",naxes[0],cptr);
      /* Status */
      fpixel[0]=1; fpixel[1]=5;
      for (i=0; i<spec[j].th.np; i++)
	dat[i]=(double)(cspec->st[spec[j].th.csidx+i]);
      if (fits_write_pix(outfits,TDOUBLE,fpixel,naxes[0],dat,&status))
	errormsg("UVES_wFITSfile(): Cannot write status array of length %d\n\
\tto FITS file %s",naxes[0],cptr);
      /** Write out keywords which relate to all spectra **/
      if (!j) {
	/* Leave mark of UVES_popler */
	sprintf(comment,
	     "Output from UVES POst Pipeline Echelle Reduction, by M. T. Murphy,");
	if (fits_write_history(outfits,comment,&status))
	  errormsg("UVES_wFITSfile(): Cannot write history keyword\n\tto file %s",
		   cspec->FITSfile);
	sprintf(comment,"UVES_popler: Version %5.2lf",VERSION);
	if (fits_write_history(outfits,comment,&status))
	  errormsg("UVES_wFITSfile(): Cannot write history keyword\n\tto file %s",
		   cspec->FITSfile);
	sprintf(comment,"\t%s",WWW);
	if (fits_write_history(outfits,comment,&status))
	  errormsg("UVES_wFITSfile(): Cannot write history keyword\n\tto file %s",
		   cspec->FITSfile);
	/* Date of writing */
	if (fits_write_date(outfits,&status))
	  errormsg("UVES_wFITSfile(): Cannot write date to file %s",cptr);
	/* Object name */
	sprintf(comment,"Object name");
	if (fits_write_key(outfits,TSTRING,"OBJECT",cspec->obj,comment,&status))
	  errormsg("UVES_wFITSfile(): Cannot write object name to file\n\t%s",cptr);
      }
      /* Heliocentric velocity of object */
      sprintf(comment,"Heliocentric velocity [km/s]");
      if (fits_write_key(outfits,TDOUBLE,"HELIOVEL",&(spec[j].vhel),comment,
			 &status))
	errormsg("UVES_wFITSfile(): Cannot write heliocentric velocity to file\n\
\t%s",cptr);
      /* Median FWHM of ThAr lines */
      sprintf(comment,"Median FWHM of ThAr lines [km/s]");
      if (fits_write_key(outfits,TDOUBLE,"ARCFWHM",&(spec[j].arcfwhm),comment,
			 &status))
	errormsg("UVES_wFITSfile(): Cannot write heliocentric velocity to file\n\
\t%s",cptr);
      /* Write out actual, original and archival ThAr filenames */
      sprintf(comment,"ThAr file name");
      if (fits_write_key(outfits,TSTRING,"THFILE",spec[j].abthfile,comment,
			 &status))
	errormsg("UVES_wFITSfile(): Could not add THFILE keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"ThAr archival file name");
      if (fits_write_key(outfits,TSTRING,"THARFILE",spec[j].tharfile,comment,
			 &status))
	errormsg("UVES_wFITSfile(): Could not add THARFILE keyword to file\n\t%s",
		 cptr);
      /* Write out actual, original and archival object filenames */
      sprintf(comment,"Object file name");
      if (fits_write_key(outfits,TSTRING,"OBFILE",spec[j].abfile,comment,&status))
	errormsg("UVES_wFITSfile(): Could not add OBFILE keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Object archival file name");
      if (fits_write_key(outfits,TSTRING,"OBARFILE",spec[j].arfile,comment,
			 &status))
	errormsg("UVES_wFITSfile(): Could not add OBARFILE keyword to file\n\t%s",
		 cptr);
      /* Write header keywords for wavelength information */
      if (par->linear) { crval1=spec[j].th.fwl; cd1_1=cspec->dwl; dc_flag=0; }
      else { crval1=log10(spec[j].th.fwl); cd1_1=cspec->dv; dc_flag=1; }
      /* Starting pixel */
      sprintf(comment,"Starting pixel (1-indexed)");
      if (fits_write_key(outfits,TINT,"CRPIX1",&crpix1,comment,&status))
	errormsg("UVES_wFITSfile(): Could not add CRPIX1 keyword to file\n\t%s",
		 cptr);
      /* Log-linear flag */
      sprintf(comment,"Log-linear flag: Linear=0, Log-linear=1");
      if (fits_write_key(outfits,TINT,"DC-FLAG",&dc_flag,comment,&status))
	errormsg("UVES_wFITSfile(): Could not add DC_FLAG keyword to file\n\t%s",
		 cptr);
      /* Starting wavelength */
      if (par->linear) sprintf(comment,"Center wavelength of first pixel");
      else sprintf(comment,"Center wavelength (log10) of first pixel");
      if (fits_write_key(outfits,TDOUBLE,"CRVAL1",&crval1,comment,&status))
	errormsg("UVES_wFITSfile(): Could not add CRVAL1 keyword to file\n\t%s",
		 cptr);
      /* Dispersion */
      if (par->linear)
	sprintf(comment,"Dispersion: w[i]=w[0]+([i+1]-CRPIX1)*CD1_1");
      else sprintf(comment,"Dispersion: w[i]=wl[0]*10^{([i+1]-CRPIX1)*CD1_1}");
      if (fits_write_key(outfits,TDOUBLE,"CD1_1",&cd1_1,comment,&status))
	errormsg("UVES_wFITSfile(): Could not add CD1_1 keyword to file\n\t%s",
		 cptr);
      if (fits_write_key(outfits,TDOUBLE,"CDELT1",&cd1_1,comment,&status))
	errormsg("UVES_wFITSfile(): Could not add CDELT1 keyword to file\n\t%s",
		 cptr);
      /* Temperature and pressure inside spectrograph */ 
      sprintf(comment,"Temperature near echelle grating for ThAr");
      if (fits_write_key(outfits,TDOUBLE,"THTEMP",&(spec->thtemp),comment,&status))
	errormsg("UVES_wFITSfile(): Could not add THTEMP keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Temperature near echelle grating for object");
      if (fits_write_key(outfits,TDOUBLE,"OBTEMP",&(spec->temp),comment,&status))
	errormsg("UVES_wFITSfile(): Could not add OBTEMP keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Pressure inside spectrograph for ThAr");
      if (fits_write_key(outfits,TDOUBLE,"THPRES",&(spec->thpres),comment,&status))
	errormsg("UVES_wFITSfile(): Could not add THPRES keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Pressure inside spectrograph for object");
      if (fits_write_key(outfits,TDOUBLE,"OBPRES",&(spec->pres),comment,&status))
	errormsg("UVES_wFITSfile(): Could not add OBPRES keyword to file\n\t%s",
		 cptr);
      /* UVES_popler-specific values */
      sprintf(comment,"Center wavelength of first pixel [Angstroms]");
      if (fits_write_key(outfits,TDOUBLE,"UP_WLSRT",&(spec[j].th.fwl),comment,
			 &status))
	errormsg("UVES_wFITSfile(): Could not add UP_WLSRT keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Center wavelength of last pixel [Angstroms]");
      if (fits_write_key(outfits,TDOUBLE,"UP_WLEND",&(spec[j]).th.lwl,comment,
			 &status))
	errormsg("UVES_wFITSfile(): Could not add UP_WLEND keyword to file\n\t%s",
		 cptr);
      if (par->linear) { disp=cd1_1; sprintf(comment,"Dispersion [Angstroms]"); }
      else { disp=par->disp; sprintf(comment,"Dispersion [km/s]"); }
      if (fits_write_key(outfits,TDOUBLE,"UP_DISP",&disp,comment,&status))
	errormsg("UVES_wFITSfile(): Could not add UP_DISP keyword to file\n\t%s",
		 cptr);
      /* Write descriptive keywords for each array */
      sprintf(comment,"Arbitrary units");
      if (fits_write_key(outfits,TSTRING,"UP_ARR01","ThAr flux spectrum",
			 comment,&status))
	errormsg("UVES_wFITSfile(): Could not add UP_ARR01 keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Arbitrary units; No info. in pixel if <0.0");
      if (fits_write_key(outfits,TSTRING,"UP_ARR02","ThAr error",comment,&status))
	errormsg("UVES_wFITSfile(): Could not add UP_ARR02 keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"1=valid; -ve integer=clipped - see web for code");
      if (fits_write_key(outfits,TSTRING,"UP_ARR03","ThAr status",comment,&status))
	errormsg("UVES_wFITSfile(): Could not add UP_ARR03 keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"As for UP_ARR02; add HELIOVEL to wavelength");
      if (fits_write_key(outfits,TSTRING,"UP_ARR04","Norm. object error",comment,
			 &status))
	errormsg("UVES_wFITSfile(): Could not add UP_ARR04 keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"1=valid; -ve integer=clipped - see web for code");
      if (fits_write_key(outfits,TSTRING,"UP_ARR05","Object status",comment,
			 &status))
	errormsg("UVES_wFITSfile(): Could not add UP_ARR05 keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Arbitrary units");
      if (fits_write_key(outfits,TSTRING,"ARRAY1","ThAr flux spectrum",
			 comment,&status))
	errormsg("UVES_wFITSfile(): Could not add ARRAY1 keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Arbitrary units; No info. in pixel if <0.0");
      if (fits_write_key(outfits,TSTRING,"ARRAY2","ThAr error",comment,&status))
	errormsg("UVES_wFITSfile(): Could not add ARRAY2 keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"1=valid; -ve integer=clipped - see web for code");
      if (fits_write_key(outfits,TSTRING,"ARRAY3","ThAr status",comment,&status))
	errormsg("UVES_wFITSfile(): Could not add ARRAY3 keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"As for ARRAY2; add HELIOVEL to wavelength");
      if (fits_write_key(outfits,TSTRING,"ARRAY4","Norm. object error",comment,
			 &status))
	errormsg("UVES_wFITSfile(): Could not add ARRAY4 keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"1=valid; -ve integer=clipped - see web for code");
      if (fits_write_key(outfits,TSTRING,"ARRAY5","Object status",comment,
			 &status))
	errormsg("UVES_wFITSfile(): Could not add ARRAY5 keyword to file\n\t%s",
		 cptr);
    }
    /* Close file */
    fits_close_file(outfits,&status);
  }

  /* Clean up */
  free(dat);

  return 1;

}
