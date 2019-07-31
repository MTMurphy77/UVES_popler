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
  double   *dtbl=NULL;
  double   version=0.0;
  long     naxes[2]={0,0};
  long     fpixel[2]={0,0};
  int      allUVES=1;
  int      crpix1=1,dc_flag=0;
  int      status=0,naxis=2,bitpix=DOUBLE_IMG;
  int      *itbl=NULL;
  int      i=0,j=0,k=0,l=0;
  char     ctype1[FLEN_KEYWORD]="\0";
  char     fitsname[LNGSTRLEN]="\0";
  char     comment[FLEN_COMMENT]="\0";
  char     key[FLEN_KEYWORD]="\0";
  char     path[LNGSTRLEN]="\0",pathend[LNGSTRLEN]="\0";
  char     delim[2]="/";
  char     *token[FLEN_KEYWORD];
  char     *extname[3]={"CombinedSpecInfo","WavelengthCoverage","ExposureInfo"};
  char *tcomt[TNCOLCOM]={"Wavelength","SNRmax","SNRmed","CNRmax","CNRmed","NomResolPower",\
			 "ArcResolPower"};
  char *tcomf[TNCOLCOM]={"1D","1D","1D","1D","1D","1D","1D"};
  char *tcomu[TNCOLCOM]={"Central wavelength of coarse information bin",\
			 "Signal-to-noise ratio maximum in bin",\
			 "Signal-to-noise ratio median in bin",\
			 "Continuum-to-noise ratio maximum in bin",\
			 "Continuum-to-noise ratio median in bin",\
			 "Median nominal resolving power in bin",\
			 "Median arc line resolving power in bin"};
  char *twcot[TNCOLCOM]={"LambdaStart","LambdaEnd"};
  char *twcof[TNCOLCOM]={"1D","1D"};
  char *twcou[TNCOLWCO]={"Starting wavelength of region (cent. of first pix)",\
			 "Ending wavelength of region (cent. of last pix)"};
  char *texpt[TNCOLEXP]={"FilePathEnd","FluxFileName","ErrFileName","WavFileName","ArchiveFileName",\
			 "WavcalArchiveFileName","FileType","ObjectName","ProgramID","OBID",\
			 "CenWav","UTDate","DPR-TECH","DPR-TYPE","DPR-CATG",\
			 "INS-PATH","INS-MODE","INS-DROT","Binning","JD",\
			 "Epoch","UTStart","ExpTime","RA","Dec",\
                         "SlitWidth","Encoder","AirmassStart","AirmassEnd","SeeingStart",\
			 "SeeingEnd","SeeingExtract","Latitude","Longitude","Altitude",\
			 "HelioVelCalc","HelioVelHead","BaryVelHead","Temperature","Pressure",\
			 "MoonRA","MoonDec","Moon_Angle","MoonPhase","WavcalUTDate",\
			 "WavcalJD","WavcalUTStart","WavcalExpTime","WavcalSlitWidth","WavcalEncoder",\
			 "WavcalTemp","WavcalPress","WavcalArcLines","WavcalRMS","WavcalArcFWHM"};
  char *texpf[TNCOLEXP]={"72A","72A","72A","72A","72A","72A","1I","64A","64A","1K",\
                         "1I","15A","64A","64A","64A","64A","64A","64A","8A","1D",\
                         "1D","1D","1D","1D","1D","1D","1K","1D","1D","1D",\
			 "1D","1D","1D","1D","1D","1D","1D","1D","1D","1D",\
			 "1D","1D","1D","1D","15A","1D","1D","1D","1D","1K",\
			 "1D","1D","1J","1D","1D"};
  char *texpu[TNCOLEXP]={" "," "," "," "," "," "," "," "," "," ",\
		       "Nominal central wavelength [nm]",\
		       "UT date at exposure start [YYY-MM-DD]",\
		       " "," "," "," "," "," ",\
		       "CCD on-chip binning (spectral x spatial)",\
		       "Julian date at exp. start [days]",\
		       "Epoch at exposure start [yrs]",\
		       "UT at exposure start [hrs]",\
		       "Exposure time [s]",\
		       "Right ascension of telescope [hrs]",\
		       "Declination of telescope [deg]",\
		       "Slit width [arcseconds]",\
		       "Grating encoder value",\
		       "Airmass at start of exposure",\
		       "Airmass at end of exposure",\
		       "FWHM seeing at start of exposure [arcseconds]",\
		       "FWHM seeing at end of exposure [arcseconds]",\
		       "FWHM seeing from extraction [arcseconds]",\
		       "Latitude of observatory [deg]",\
		       "Longitude: positive=West [deg]",\
		       "Altitude of observatory [m]",\
		       "Heliocentric velocity (computed) [km/s]",\
		       "Heliocentric velocity (header) [km/s]",\
		       "Barycentric velocity (header) [km/s]",\
		       "Temperature at grating [deg C or K]",\
		       "Av. atmospheric pressure [hPa]",\
		       "Moon RA [hrs]",\
		       "Moon Dec. [deg]",\
		       "Angle between object and moon [deg]",\
		       "Moon phase as fraction of period",\
		       "UT date of wavcal start [YYY-MM-DD]",\
		       "Julian date for wavcal [days]",\
		       "UT at wavcal exposure start [hrs]",\
                       "Wavcal exposure time [s]",\
		       "Slit width for wavcal [arcseconds]",\
		       "Grating encoder value for wavcal",\
		       "Temp. at grating for ThAr [deg C or K]",\
		       "Av. atmosph. pressure for ThAr [hPa]",\
		       "Number of arc lines used in wavcal.",\
		       "RMS of wavcal [m/s]",\
		       "Median FWHM of ThAr lines [km/s]"};
  char     **stbl;
  char     *cptr;
  fitsfile *outfits=NULL;

#define ERR_PRIKEY { \
    errormsg("UVES_wFITSfile(): Could not add\n\
\t'%s' keyword to primary header of FITS file\n\t%s",key,cspec->FITSfile); }
#define ERR_TABCOM { \
    errormsg("UVES_wFITSfile(): Cannot add column\n\
\t'%s' to binary table extension %d, '%s', in FITS file\n\t%s",\
	     tcomt[j],0,extname[0],cspec->FITSfile); }
#define ERR_TABWCO { \
    errormsg("UVES_wFITSfile(): Cannot add column\n\
\t'%s' to binary table extension %d, '%s', in FITS file\n\t%s",\
	     twcot[j],1,extname[1],cspec->FITSfile); }
#define ERR_TABEXP { \
    errormsg("UVES_wFITSfile(): Cannot add column\n\
\t'%s' to binary table extension %d, '%s', in FITS file\n\t%s",\
	     texpt[j],2,extname[2],cspec->FITSfile); }

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
  /* Date of writing */
  if (fits_write_date(outfits,&status))
    errormsg("UVES_wFITSfile(): Cannot write date to file %s",cspec->FITSfile);
  /* Object name */
  sprintf(key,"OBJECT"); sprintf(comment,"Object name");
  if (fits_write_key(outfits,TSTRING,key,cspec->obj,comment,&status)) { ERR_PRIKEY; }
  /* Write header keywords for wavelength information */
  if (par->linear) {
    crval1=cspec->wl[0]; cd1_1=cspec->dwl; dc_flag=0; sprintf(ctype1,"LINEAR");
  } else {
    crval1=log10(cspec->wl[0]); cd1_1=cspec->dv; dc_flag=1; sprintf(ctype1,"LOGLIN");
  }
  /* Starting pixel */
  sprintf(key,"CRPIX1"); sprintf(comment,"Starting pixel (1-indexed)");
  if (fits_write_key(outfits,TINT,key,&crpix1,comment,&status)) { ERR_PRIKEY; }
  /* Log-linear flag */
  sprintf(key,"CTYPE1");
  sprintf(comment,"Log-lin flag: Lin='LINEAR', Log-lin='LOGLIN'");
  if (fits_write_key(outfits,TSTRING,key,ctype1,comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"DC-FLAG");
  sprintf(comment,"Log-linear flag: Linear=0, Log-linear=1");
  if (fits_write_key(outfits,TINT,key,&dc_flag,comment,&status)) { ERR_PRIKEY; }
  /* Starting wavelength */
  sprintf(key,"CRVAL1");
  if (par->linear) sprintf(comment,"Center wavelength of first pixel");
  else sprintf(comment,"Center wavelength (log10) of first pixel");
  if (fits_write_key(outfits,TDOUBLE,key,&crval1,comment,&status)) { ERR_PRIKEY; }
  /* Dispersion */
  sprintf(key,"CD1_1");
  if (par->linear) sprintf(comment,"Dispersion: w[i]=w[0]+([i+1]-CRPIX1)*CD1_1");
  else sprintf(comment,"Dispersion: w[i]=wl[0]*10^{([i+1]-CRPIX1)*CD1_1}");
  if (fits_write_key(outfits,TDOUBLE,key,&cd1_1,comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"CDELT1");
  if (fits_write_key(outfits,TDOUBLE,key,&cd1_1,comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"CD2_2"); sprintf(comment,"");
  if (fits_write_key(outfits,TDOUBLE,key,&cd2_2,comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"CDELT2");
  if (fits_write_key(outfits,TDOUBLE,key,&cd2_2,comment,&status)) { ERR_PRIKEY; }
  /* UVES_popler-specific values */
  sprintf(key,"UP_WLSRT");
  sprintf(comment,"Center wavelength of first pixel [Angstroms]");
  if (fits_write_key(outfits,TDOUBLE,key,&cspec->wl[0],comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"UP_WLEND");
  sprintf(comment,"Center wavelength of last pixel [Angstroms]");
  if (fits_write_key(outfits,TDOUBLE,key,&cspec->wl[cspec->np-1],comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"UP_ZEM"); sprintf(comment,"Emission redshift");
  if (fits_write_key(outfits,TDOUBLE,key,&par->zem,comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"UP_DISP");
  if (par->linear) { disp=cd1_1; sprintf(comment,"Dispersion [Angstroms]"); }
  else { disp=par->disp; sprintf(comment,"Dispersion [km/s]"); }
  if (fits_write_key(outfits,TDOUBLE,key,&disp,comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"UP_NEXP"); sprintf(comment,"Number of exposures");
  if (fits_write_key(outfits,TINT,key,&cspec->nexp,comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"UP_TEXP"); sprintf(comment,"Total exposure time [s]");
  if (fits_write_key(outfits,TDOUBLE,key,&cspec->texp,comment,&status))
    { ERR_PRIKEY; }
  /* Write descriptive keywords for each array */
  sprintf(key,"UP_ARR01"); sprintf(comment,"No units");
  if (fits_write_key(outfits,TSTRING,key,"Normalized flux spectrum",
		     comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"UP_ARR02"); sprintf(comment,"No units; No info. in pixel if <0.0");
  if (fits_write_key(outfits,TSTRING,key,"Normalized error",
		     comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"UP_ARR03");
  if (fits_write_key(outfits,TSTRING,key,"Norm. expected fluctuation",
		     comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"UP_ARR04"); sprintf(comment,"Arbitrary units");
  if (fits_write_key(outfits,TSTRING,key,"Continuum by which data are normalized",
		     comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"UP_ARR05");
  sprintf(comment,"1=valid; -ve integer=clipped - see web for code");
  if (fits_write_key(outfits,TSTRING,key,"Status",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"UP_ARR06");
  sprintf(comment,"# pix contrib. to comb. pix. before sigma-clip");
  if (fits_write_key(outfits,TSTRING,key,"NPIX before clip",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"UP_ARR07");
  sprintf(comment,"# pix contrib. to comb. pix. after sigma-clip");
  if (fits_write_key(outfits,TSTRING,key,"NPIX after clip",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"UP_ARR08");
  sprintf(comment,"chisq. about weighted mean before sigma-clip");
  if (fits_write_key(outfits,TSTRING,key,"chisq. before clip",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"UP_ARR09");
  sprintf(comment,"chisq. about weighted mean after sigma-clip");
  if (fits_write_key(outfits,TSTRING,key,"chisq. after clip",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"ARRAY1"); sprintf(comment,"No units");
  if (fits_write_key(outfits,TSTRING,key,"Normalized flux spectrum",comment,
		     &status)) { ERR_PRIKEY; }
  sprintf(key,"ARRAY2"); sprintf(comment,"No units; No info. in pixel if <0.0");
  if (fits_write_key(outfits,TSTRING,key,"Normalized error",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"ARRAY3");
  if (fits_write_key(outfits,TSTRING,key,"Norm. expected fluctuation",
		     comment,&status)) { ERR_PRIKEY; }
  sprintf(key,"ARRAY4"); sprintf(comment,"Arbitrary units");
  if (fits_write_key(outfits,TSTRING,key,
		     "Continuum by which data are normalized",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"ARRAY5");
  sprintf(comment,"1=valid; -ve integer=clipped - see web for code");
  if (fits_write_key(outfits,TSTRING,key,"Status",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"ARRAY6");
  sprintf(comment,"# pix contrib. to comb. pix. before sigma-clip");
  if (fits_write_key(outfits,TSTRING,key,"NPIX before clip",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"ARRAY7");
  sprintf(comment,"# pix contrib. to comb. pix. after sigma-clip");
  if (fits_write_key(outfits,TSTRING,key,"NPIX after clip",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"ARRAY8");
  sprintf(comment,"chisq. about weighted mean before sigma-clip");
  if (fits_write_key(outfits,TSTRING,key,"chisq. before clip",comment,&status))
    { ERR_PRIKEY; }
  sprintf(key,"ARRAY9");
  sprintf(comment,"chisq. about weighted mean after sigma-clip");
  if (fits_write_key(outfits,TSTRING,key,"chisq. after clip",comment,&status))
    { ERR_PRIKEY; }

  /* Write out additional information about the combined and
     contributing exposures in binary table extensions, but only do
     this for UVES spectra */
  for (i=0; i<nspec; i++) if (spec[i].ftype!=FTUVES) { allUVES=0; break; }
  if (allUVES) {
    /** Create binary table HDU containing information about the combined spectrum **/
    /* Create table */
    if (fits_create_tbl(outfits,BINARY_TBL,0,TNCOLCOM,tcomt,tcomf,tcomu,extname[0],&status))
      errormsg("UVES_wFITSfile(): Cannot create binary table extension\n\
\t'%s' in FITS file\n\t%s",extname[0],cspec->FITSfile);
    /* Write coarse information arrays */
    j=0;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,cspec->nc,cspec->cwav,&status)) { ERR_TABCOM; }
    j++;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,cspec->nc,cspec->csnrmax,&status)) { ERR_TABCOM; }
    j++;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,cspec->nc,cspec->csnrmed,&status)) { ERR_TABCOM; }
    j++;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,cspec->nc,cspec->ccnrmax,&status)) { ERR_TABCOM; }
    j++;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,cspec->nc,cspec->ccnrmed,&status)) { ERR_TABCOM; }
    j++;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,cspec->nc,cspec->cresnom,&status)) { ERR_TABCOM; }
    j++;
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,cspec->nc,cspec->cresarc,&status)) { ERR_TABCOM; }

    /** Create binary table HDU containing the wavelength coverage map **/
    /* Create table */
    if (fits_create_tbl(outfits,BINARY_TBL,0,TNCOLWCO,twcot,twcof,twcou,extname[1],&status))
      errormsg("UVES_wFITSfile(): Cannot create binary table extension\n\
\t'%s' in FITS file\n\t%s",extname[1],cspec->FITSfile);
    /* Allocate memory for temporary double array */
    if ((dtbl=darray(cspec->nwavcov))==NULL)
      errormsg("UVES_wFITSfile(): Cannot allocate memory for dtbl\n\
\tarray of length %d",cspec->nwavcov);
    /* Write wavelength coverage arrays */
    j=0; for (i=0; i<cspec->nwavcov; i++) dtbl[i]=cspec->wavcov[i][0];
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,cspec->nwavcov,dtbl,&status)) { ERR_TABWCO; }
    j++; for (i=0; i<cspec->nwavcov; i++) dtbl[i]=cspec->wavcov[i][1];
    if (fits_write_col(outfits,TDOUBLE,j+1,1,1,cspec->nwavcov,dtbl,&status)) { ERR_TABWCO; }

    /** Create binary table HDU containing parameters from each exposure **/
    /* Only do this if we have contributing spectra, i.e. not for
       already-combined spectra */
    if (nspec) {
      /* Create table */
      if (fits_create_tbl(outfits,BINARY_TBL,0,TNCOLEXP,texpt,texpf,texpu,
			  extname[2],&status))
	errormsg("UVES_wFITSfile(): Cannot create binary table extension\n\
\t'%s' in FITS file\n\t%s",extname[2],cspec->FITSfile);
      /* Write number of files to header (may not be number of exposures) */
      sprintf(key,"NFILES"); sprintf(comment,"Number of flux files");
      if (fits_write_key(outfits,TINT,key,&nspec,comment,&status)) {
	errormsg("UVES_wFITSfile(): Could not add\n\
\t'%s' keyword to binary table extension %d, '%s', in FITS file\n\t%s",
		 key,2,extname[2],cspec->FITSfile);
      }
      /* Allocate memory for double and integer arrays */
      if ((dtbl=darray(nspec))==NULL)
	errormsg("UVES_wFITSfile(): Cannot allocate memory for dtbl\n\
\tarray of length %d",nspec);
      if ((itbl=iarray(nspec))==NULL)
	errormsg("UVES_wFITSfile(): Cannot allocate memory for itbl\n\
\tarray of length %d",nspec);
      /* Allocate memory for string matrix */
      if ((stbl=cmatrix(nspec,FLEN_KEYWORD))==NULL)
	errormsg("UVES_wFITSfile(): Cannot allocate memory for stbl\n\
\tarray of length %dx%d",nspec,FLEN_KEYWORD);
      /* Paths: Last 3 elements of path only */
      j=0;
      for (i=0; i<nspec; i++) {
	/* Determine last 3 elements of path */
	sprintf(path,"%s",spec[i].path); token[0]=strtok(path,delim);
	k=0; while (token[k]!=NULL) { k++; token[k]=strtok(NULL,delim); }
	sprintf(pathend,"/"); for (l=(MAX(0,k-3)); l<k; l++) {
	  strcat(pathend,token[l]); strcat(pathend,"/");
	}
	if (strlen(pathend)<FLEN_KEYWORD) sprintf(stbl[i],"%s",pathend);
	else errormsg("UVES_wFITSfile(): Length of last 3 elements of file path,\n\
\t%s,\n\too long (=%d) to be written to output FITS file\n\t%s.\n\
\tConsider decreasing path name length and/or increasing length of\n\
\tinternal variable beyond FLEN_KEYWORD (=%d)",pathend,strlen(pathend),
		      cspec->FITSfile,FLEN_KEYWORD);
      }
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      /* Flux, error and wavelength file names */
      j++;
      for (i=0; i<nspec; i++) {
	if (strlen(spec[i].abfile)<FLEN_KEYWORD) sprintf(stbl[i],"%s",spec[i].abfile);
	else errormsg("UVES_wFITSfile(): Length of flux file name,\n\t%s,\n\
\ttoo long (=%d) to be written to output FITS file\n\t%s.\n\
\tConsider decreasing file name length and/or increasing length of\n\
\tinternal variable beyond FLEN_KEYWORD (=%d)",spec[i].abfile,
		      strlen(spec[i].abfile),cspec->FITSfile,FLEN_KEYWORD);
      }
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      for (i=0; i<nspec; i++) {
	if (strlen(spec[i].aberfile)<FLEN_KEYWORD) sprintf(stbl[i],"%s",spec[i].aberfile);
	else errormsg("UVES_wFITSfile(): Length of error file name,\n\t%s,\n\
\ttoo long (=%d) to be written to output FITS file\n\t%s.\n\
\tConsider decreasing file name length and/or increasing length of\n\
\tinternal variable beyond FLEN_KEYWORD (=%d)",spec[i].aberfile,
		      strlen(spec[i].aberfile),cspec->FITSfile,FLEN_KEYWORD);
      }
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      for (i=0; i<nspec; i++) {
	if (strlen(spec[i].abthfile)<FLEN_KEYWORD) sprintf(stbl[i],"%s",spec[i].abwlfile);
	else errormsg("UVES_wFITSfile(): Length of wavelength file name,\n\t%s,\n\
\ttoo long (=%d) to be written to output FITS file\n\t%s.\n\
\tConsider decreasing file name length and/or increasing length of\n\
\tinternal variable beyond FLEN_KEYWORD (=%d)",spec[i].abwlfile,
		      strlen(spec[i].abwlfile),cspec->FITSfile,FLEN_KEYWORD);
      }
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* Archival file name */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].arfile);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal archival file name */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].wlarfile);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* File type */
      for (i=0; i<nspec; i++) itbl[i]=spec[i].ftype;
      if (fits_write_col(outfits,TINT,j+1,1,1,nspec,itbl,&status)) { ERR_TABEXP; }
      j++;
      /* Object name */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].obj);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* Program ID */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].progid);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* OB ID */
      for (i=0; i<nspec; i++) itbl[i]=spec[i].obid;
      if (fits_write_col(outfits,TINT,j+1,1,1,nspec,itbl,&status)) { ERR_TABEXP; }
      j++;
      /* Nominal central wavelength */
      for (i=0; i<nspec; i++) itbl[i]=rint(spec[i].cwl);
      if (fits_write_col(outfits,TINT,j+1,1,1,nspec,itbl,&status)) { ERR_TABEXP; }
      j++;
      /* Date */
      for (i=0; i<nspec; i++)
	sprintf(stbl[i],"%04d-%02d-%02d",spec[i].year,spec[i].month,spec[i].day);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* DPR.TECH */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].dprtech);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* DPR.TYPE */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].dprtype);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* DPR.CATG */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].dprcatg);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* INS.PATH */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].inspath);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* INS.MODE */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].insmode);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* INS.DROT */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%s",spec[i].insdrot);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* Binning */
      for (i=0; i<nspec; i++) sprintf(stbl[i],"%dx%d",spec[i].binx,spec[i].biny);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* Julian date */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].jd;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Epoch */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].epoch;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* UT */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].ut;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Exposure time */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].etime;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* RA */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].ra;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* DEC */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].dec;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Slit width */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].sw;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Grating encoder value */
      for (i=0; i<nspec; i++) itbl[i]=spec[i].encoder;
      if (fits_write_col(outfits,TINT,j+1,1,1,nspec,itbl,&status)) { ERR_TABEXP; }
      j++;
      /* Airmass */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].airmass[0];
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].airmass[1];
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Seeing */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].seeing[0];
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].seeing[1];
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].seeing[2];
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Latitude */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].lat;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Longitude */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].lon;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Altitude */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].alt;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Heliocentric velocity (calculated) */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].vhel;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Heliocentric velocity (from original file header) */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].vhel_head;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Barycentric velocity (from original file header) */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].vbar_head;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Temperature */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].temp;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Pressure */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].pres;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Moon RA */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].moon_ra;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Moon Dec */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].moon_dec;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Moon angular distance */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].moon_ang;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Moon phase */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].moon_phase;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal date */
      for (i=0; i<nspec; i++)
	sprintf(stbl[i],"%04d-%02d-%02d",spec[i].wc_year,spec[i].wc_month,spec[i].wc_day);
      if (fits_write_col(outfits,TSTRING,j+1,1,1,nspec,stbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal Julian date */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].wc_jd;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal UT */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].wc_ut;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal exposure time */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].wc_etime;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal slit width */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].wc_sw;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal grating encoder value */
      for (i=0; i<nspec; i++) itbl[i]=spec[i].wc_encoder;
      if (fits_write_col(outfits,TINT,j+1,1,1,nspec,itbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal temperature */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].wc_temp;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal pressure */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].wc_pres;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal number of lines in wavelength solution */
      for (i=0; i<nspec; i++) itbl[i]=spec[i].ts.np;
      if (fits_write_col(outfits,TINT,j+1,1,1,nspec,itbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal median RMS of wavelength solution */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].wc_resid;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Wavcal arc FWHM */
      for (i=0; i<nspec; i++) dtbl[i]=spec[i].arcfwhm;
      if (fits_write_col(outfits,TDOUBLE,j+1,1,1,nspec,dtbl,&status)) { ERR_TABEXP; }
      j++;
      /* Clean up */
      free(dtbl); free(itbl); free(*stbl); free(stbl);
    }
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
      if (fits_write_key(outfits,TDOUBLE,"THTEMP",&(spec->wc_temp),comment,&status))
	errormsg("UVES_wFITSfile(): Could not add THTEMP keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Temperature near echelle grating for object");
      if (fits_write_key(outfits,TDOUBLE,"OBTEMP",&(spec->temp),comment,&status))
	errormsg("UVES_wFITSfile(): Could not add OBTEMP keyword to file\n\t%s",
		 cptr);
      sprintf(comment,"Pressure inside spectrograph for ThAr");
      if (fits_write_key(outfits,TDOUBLE,"THPRES",&(spec->wc_pres),comment,&status))
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

/* 
Other information to store?:

For each fxb_* file read in:
- number of echelle orders, stored as spec.nor
- number of pixels per order, stored as spec.or[0].np

For each order?:
- Order number
- Diffraction order number
- Start and end wavelengths, before cuts applied
- Final start and end wavelengths, after cuts applied
- Seeing derived from optimal extraction
- Number of ThAr lines used in solution
- RMS residuals of ThAr lines used in solution
- Median FWHM resolution from wavelength solution
- Median SNR and CNR near middle of order

*/
