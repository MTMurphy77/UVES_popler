/****************************************************************************

UVES_popler: UVES POst PIpeline Echelle Reduction

****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "UVES_popler.h"
#include "stats.h"
#include "input.h"
#include "error.h"

/* Global declarations */
char      *progname;
plotenv   plenv;

/****************************************************************************
* Print the usage message
****************************************************************************/

void usage(void) {

  fprintf(stderr,"%s (UVES POst PipeLine Echelle Reduction) by Michael Murphy\n\
 %s\n Version: %.2lf (%s)\n\n",progname,WWW,VERSION,DATECREATE);

  fprintf(stderr,"Usage: %s [OPTIONS] [FLUX FILE LIST or UPL file]\n\n", progname);
  fprintf(stderr, "Options:\n\
General spectrum specifications:\n\
 -combmeth   =       %1d    : Comb. method, order-scaling (0), order-fitting (1)\n\
 -disp       =       TBD  : Dispersion in km/s or Angstroms (for lin. disp.)\n\
 -filetype   =       %1d    : File origin: UVES=0, IRAF=1, MAKEE=2, IRAFLSS=3,\n\
                               HIREDUX=4, ESOMERGED=5, MAGE=6, IRAFESI=7,\n\
                               COMB=9, mixed=-1\n\
 -helio      =       %1d    : Input files in observed(0) or heliocentric(1) frame\n\
 -vacwl      =       %1d    : Input files in air (0) or vacuum (1) wavelengths\n\
 -zem        = %10.5lf : Emission redshift of QSO (if required for cont. fit)\n\
 -linear     =       %1d    : Linear (1) or log-linear (0) dispersion\n\
 -thar       =       %1d    : ThAr spectra: Don't include (0), re-bin and write\n\
                               out (1) or combine (2)\n\
\n\
Initial sigma-clipping of echelle orders (not done when using -thar):\n\
 -nordclip   =      %2d    : No. most positive pixels removed from clip window\n\
 -nordsig    =      %2d    : No. pixels in order sig.-clip window\n\
 -ordsig     = %10.2lf : Sigma-clipping level for orders [rms in window]\n\
 -ordsignbr  = %10.2lf : Sigma-clipping level for neighbouring pixels\n\
 -ordsigzero = %10.2lf : Only clip when mean is this many sigma above zero\n\
 -ordmedfrac = %10.2lf : Median filter width relative to order length\n\
 -ordmedrej  = %10.2lf : Reject pixels this much smaller than median\n\
\n\
Spectrum/order combination parameters:\n\
 -clipsig    = %10.2lf : Sigma-clipping level for spectra combination\n\
 -lrgerr     = %10.2lf : Large error rejection ratio for combination\n\
 -rankspec   =       %1d    : (0) Autoscale; (1) Hi-SNR ord of 1st spec unscaled\n\
 -scalmeth   =       %1d    : Ord. scal. meth.: (0) Ratio+sig.clip, (1) Chisq.min.\n\
 -scalclip   = %10.2lf : Sigma-clipping level for inter-order scale-factor\n\
 -scalerr    = %10.2lf : Limit for relative error on inter-order scale-factor\n\
 -nscalclip  =       %1d    : Max. rejection iterations for scalmeth=1\n\
\n\
Order continuum fitting method parameters:\n\
 -contftyp   =      %2d    : Type of fit (1=polynomial, 2=Chebyshev, 3=Legendre)\n\
 -contord    =      %2d    : Order of fitted continuum (1=const.)\n\
 -contpctl   = %10.2lf : Lower percentile to remove for initial fit [%%/100]\n\
 -contsigl   = %10.2lf : Lower rej. level for order continuum fitting [sigma]\n\
 -contsigu   = %10.2lf : Upper rej. level for order continuum fitting [sigma]\n\
 -contwgt    =    %4d    : Weighting scale in pixels for ends of orders\n\
\n\
Combined spectrum continuum fitting method parameters:\n\
 -cordlya    =      %2d    : Order of fitted continuum below Lyman-alpha\n\
 -cordred    =      %2d    : Order of fitted continuum above Lyman-alpha\n\
 -ftyplya    =      %2d    : Type of fit below Lyman-alpha (1=pol, 2=Che, 3=Leg)\n\
 -ftypred    =      %2d    : Type of fit above Lyman-alpha (1=pol, 2=Che, 3=Leg)\n\
 -nocont     =       %1d    : Do not attempt automatic continuum fit\n\
 -pctllya    = %10.2lf : Lower percentile to remove for first Lya fit [%%/100]\n\
 -pctlred    = %10.2lf : Lower percentile to remove for first red fit [%%/100]\n\
 -rsiglyal   = %10.2lf : Lower rej. level for Lya continuum fit [sigma]\n\
 -rsiglyau   = %10.2lf : Upper rej. level for Lya continuum fit [sigma]\n\
 -rsigredl   = %10.2lf : Lower rej. level for cont. fit to red of Lya [sigma]\n\
 -rsigredu   = %10.2lf : Upper rej. level for cont. fit to red of Lya [sigma]\n\
 -vclya      = %10.2lf : Velocity scale for continuum fit below Lya [km/s]\n\
 -vcred      = %10.2lf : Velocity scale for continuum fit above Lya [km/s]\n\
 -vlya       = %10.2lf : Velocity below Lya where red cont. fit begins [km/s]\n\
\n\
Data file options:\n\
 -dat        =       %1d    : Write out data file with normalized spectrum\n\
 -raw        =       %1d    : Write out raw data for all orders to a FITS file\n\
 -save       =       %1d    : (1) Force output file save upon user quit;\n\
                            (2) Save+quit without display; UPL file's save flag=0)\n\
 -prefix [PREFIX]         : Specify prefix for output files, e.g. [PREFIX].upl etc.\n\
\n\
Mode options:\n\
 -replay     =       %1d    : Use interactive action replay mode\n\
 -macmap [FILE]           : Read file to map pipeline product path and file names\n\
                             between case-sensitive and case-insensitive operating\n\
                             systems, e.g. Mac and linux.\n\
 -vshift [FILE]           : Read file to apply a velocity shift(+distortion) shift\n\
                             to each spectrum. File specifies 4 cols: file path,\n\
                             velocity shift [km/s], slope [km/s/1000Ang], reference\n\
                             wavelength [Ang].\n\
 -atmomask [FILE]         : Read file to mask spectra for atmospheric/telluric\n\
                             features. File should have format: wav. start, wav. end\n\
                             [Ang.], residual intensity, wav. center [Ang.].\n\
                             Wavelengths are in vacuum.\n\
 -distort                 : Apply a random velocity shift, drift and distortion to\n\
                             each echelle order. Use this flag to remove previously\n\
                             applied distortions saved in UPL file.\n\
 -h, -help                : Print this message.\n\n",
	  COMBMETH,FILETYPE,HELIO,VACWL,ZEM,LINEAR,THAR,NORDCLIP,NORDSIG,ORDSIG,
	  ORDSIGNBR,ORDSIGZERO,ORDMEDFRAC,ORDMEDREJ,CLIPSIG,LRGERR,RANKSPEC,SCALMETH,
	  SCALCLIP,SCALERR,NSCALCLIP,CONTFTYP,CONTORD,CONTPCTL,CONTSIGL,CONTSIGU,
	  CONTWGT,CORDLYA,CORDRED,FTYPLYA,FTYPRED,NOCONT,PCTLLYA,PCTLRED,RSIGLYAL,
	  RSIGLYAU,RSIGREDL,RSIGREDU,VCLYA,VCRED,VLYA,DAT,RAW,SAVE,REPLAY);
  exit(3);
}

/****************************************************************************
* The main program

  Things to do:

****************************************************************************/

int main(int argc, char *argv[]) {
  
  int       nspec=0;     /* Number of 2D spectra read in from PWD */
  int       nact=0;      /* Number of historical actions */
  int       nact_save=0,cact=0,status=0;
  int       i=0,j=0;
  char      infile[LNGSTRLEN]="\0";
  char      query[QUERYLEN]="\0";
  char      *cptr;
  action    *act=NULL;   /* Array of historical actions to be carried out */
  atmask    amsk;        /* Atmospheric feature mask */
  cspectrum cspec;       /* Combined spectrum */
  cplot     cp;          /* Combined plot parameters */
  rplot     rp;          /* Replay plot parameters */
  params    par;         /* Parameters for reduction */
  spectrum  *spec=NULL;  /* Array of spectra */

  /* Define the program name from the command line input */
  progname=((progname=strrchr(argv[0],'/'))==NULL) ? argv[0] : progname+1;
  /* Must be at least one argument */
  if (argc==1) usage();
  /* Initialize reduction parameters */
  if (!UVES_params_init(&par))
    errormsg("Unknown error returned from UVES_params_init()");
  /* Scan command line for options */
  while (++i<argc) {
    if (!strcmp(argv[i],"-help") || !strcmp(argv[i],"-h")) usage();
    else if (!strcmp(argv[i],"-combmeth")) par.combmeth=1;
    else if (!strcmp(argv[i],"-disp")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.disp)!=1) usage();
    }
    else if (!strcmp(argv[i],"-filetype")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.filetype)!=1) usage();
    }
    else if (!strcmp(argv[i],"-helio")) par.helio=1;
    else if (!strcmp(argv[i],"-vacwl")) par.vacwl=1;
    else if (!strcmp(argv[i],"-zem")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.zem)!=1) usage();
    }
    else if (!strcmp(argv[i],"-linear")) par.linear=1;
    else if (!strcmp(argv[i],"-thar")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.thar)!=1) usage();
    }
    else if (!strcmp(argv[i],"-ordsig")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.ordsig)!=1) usage();
    }
    else if (!strcmp(argv[i],"-ordsignbr")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.ordsignbr)!=1) usage();
    }
    else if (!strcmp(argv[i],"-ordsigzero")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.ordsigzero)!=1) usage();
    }
    else if (!strcmp(argv[i],"-ordmedfrac")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.ordmedfrac)!=1) usage();
    }
    else if (!strcmp(argv[i],"-ordmedrej")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.ordmedrej)!=1) usage();
    }
    else if (!strcmp(argv[i],"-nordclip")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.nordclip)!=1) usage();
    }
    else if (!strcmp(argv[i],"-nordsig")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.nordsig)!=1) usage();
    }
    else if (!strcmp(argv[i],"-clipsig")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.clipsig)!=1) usage();
    }
    else if (!strcmp(argv[i],"-lrgerr")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.lrgerr)!=1) usage();
    }
    else if (!strcmp(argv[i],"-scalmeth")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.scalmeth)!=1) usage();
    }
    else if (!strcmp(argv[i],"-scalclip")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.scalclip)!=1) usage();
    }
    else if (!strcmp(argv[i],"-scalerr")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.scalerr)!=1) usage();
    }
    else if (!strcmp(argv[i],"-nscalclip")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.nscalclip)!=1) usage();
    }
    else if (!strcmp(argv[i],"-contpctl")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.contpctl)!=1) usage();
    }
    else if (!strcmp(argv[i],"-contsigl")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.contsigl)!=1) usage();
    }
    else if (!strcmp(argv[i],"-contsigu")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.contsigu)!=1) usage();
    }
    else if (!strcmp(argv[i],"-contftyp")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.contftyp)!=1) usage();
    }
    else if (!strcmp(argv[i],"-contord")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.contord)!=1) usage();
    }
    else if (!strcmp(argv[i],"-contwgt")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.contwgt)!=1) usage();
    }
    else if (!strcmp(argv[i],"-cordlya")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.cordlya)!=1) usage();
    }
    else if (!strcmp(argv[i],"-cordred")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.cordred)!=1) usage();
    }
    else if (!strcmp(argv[i],"-ftyplya")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.ftyplya)!=1) usage();
    }
    else if (!strcmp(argv[i],"-ftypred")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.ftypred)!=1) usage();
    }
    else if (!strcmp(argv[i],"-nocont")) par.nocont=CNTNON;
    else if (!strcmp(argv[i],"-pctllya")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.pctllya)!=1) usage();
    }
    else if (!strcmp(argv[i],"-pctlred")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.pctlred)!=1) usage();
    }
    else if (!strcmp(argv[i],"-rankspec")) par.rankspec=1;
    else if (!strcmp(argv[i],"-rsiglyal")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.rsiglyal)!=1) usage();
    }
    else if (!strcmp(argv[i],"-rsiglyau")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.rsiglyau)!=1) usage();
    }
    else if (!strcmp(argv[i],"-rsigredl")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.rsigredl)!=1) usage();
    }
    else if (!strcmp(argv[i],"-rsigredu")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.rsigredu)!=1) usage();
    }
    else if (!strcmp(argv[i],"-vclya")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.vclya)!=1) usage();
    }
    else if (!strcmp(argv[i],"-vcred")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.vcred)!=1) usage();
    }
    else if (!strcmp(argv[i],"-vlya")) {
      if (++i==argc || sscanf(argv[i],"%lf",&par.vlya)!=1) usage();
    }
    else if (!strcmp(argv[i],"-dat")) par.dat=1;
    else if (!strcmp(argv[i],"-raw")) par.raw=1;
    else if (!strcmp(argv[i],"-save")) {
      if (++i==argc || sscanf(argv[i],"%d",&par.save)!=1) usage();
    }
    else if (!strcmp(argv[i],"-prefix")) {
      if (++i==argc || !strncmp(argv[i],"-",1)) usage();
      if (sscanf(argv[i],"%s",par.prefix)!=1) usage();
    }
    else if (!strcmp(argv[i],"-replay")) par.replay=1;
    else if (!strcmp(argv[i],"-macmap")) {
      par.macmap=1; if (++i==argc || !strncmp(argv[i],"-",1)) usage();
      if (sscanf(argv[i],"%s",par.macmapfile)!=1) usage();
    }
    else if (!strcmp(argv[i],"-vshift")) {
      par.vshift=1; if (++i==argc || !strncmp(argv[i],"-",1)) usage();
      if (sscanf(argv[i],"%s",par.vshiftfile)!=1) usage();
    }
    else if (!strcmp(argv[i],"-atmomask")) {
      par.atmask=1; if (++i==argc || !strncmp(argv[i],"-",1)) usage();
      if (sscanf(argv[i],"%s",par.atmaskfile)!=1) usage();
    }
    else if (!strcmp(argv[i],"-distort")) par.distort=1;
    else if (!access(argv[i],R_OK)) {
      if (strlen(argv[i])<=LNGSTRLEN) strcpy(infile,argv[i]);
      else errormsg("Input file name too long: %s",argv[i]);
    }
    else errormsg("File %s does not exist",argv[i]);
  }
  
  /* Make sure an input file was specified */
  if (!strncmp(infile,"\0",1)) usage();

  /* Read input file */
  if (!UVES_rinputfile(infile,&spec,&nspec,&act,&nact,&cspec,&par))
    errormsg("Unknown error returned from UVES_rinputfile()");

  /* Set any unset reduction parameters */
  if (!UVES_params_set(&par))
    errormsg("Unknown error returned from UVES_params_set()");
  /* Make sure nordsig is odd */
  if (!isodd(par.nordsig)) errormsg("Nordsig (=%d) must be odd",par.nordsig);
  /* Make sure nordsig and nordclip are consistent */
  if (par.nordsig<par.nordclip+2)
    errormsg("nordsig (=%d) must be at least 2 greater\n\
\tthan nordclip (=%d)",par.nordsig,par.nordclip);
  /* Make sure nordsig is less than minimum useful order length */
  if (par.nordsig>=MINUSE)
    errormsg("nordsig (=%d) must be less than minimum useful order\n\
\tlength (=%d). Try increasing MINUSE in UVES_popler.h",par.nordsig,MINUSE);
  /* Make sure save flag is meaningful */
  if (par.save<0 || par.save>2)
    errormsg("The save flag (=%d) must be 0, 1 or 2",par.save);

  /* Read velocity shift file provided by the user */
  if (par.vshift==1) {
    if (!UVES_rvshift(&spec,nspec,&par))
      errormsg("Error reading velocity shift file %s",par.vshiftfile);
  }

  /* Read the atmospheric feature mask file provided by the user */
  if (par.atmask==1) {
    if (!UVES_ratmask(&amsk,&par))
      errormsg("Error reading atmospheric mask file %s",par.atmaskfile);
  }

  /* Read in spectra */
  for (i=0; i<nspec; i++) {
    switch (spec[i].ftype) {
    case FTUVES:
      if (!UVES_r2Dspec(&(spec[i]),&par))
	errormsg("Error returned from UVES_r2Dspec()\n\
\twhen attempting to read in spectrum %d,\n\t%s",i+1,spec[i].file);
      break;
    case FTIRAF:
      if (!UVES_r2Dspec_iraf(&(spec[i]),&par))
	errormsg("Unknown error returned from UVES_r2Dspec_iraf()\n\
\twhen attempting to read in spectrum %d,\n\t%s",i+1,spec[i].file);
      break;
    case FTMAKE:
      if (!UVES_r2Dspec_makee(&(spec[i]),&par))
	errormsg("Unknown error returned from UVES_r2Dspec_makee()\n\
\twhen attempting to read in spectrum %d,\n\t%s",i+1,spec[i].file);
      break;
    case FTIRLS:
      if (!UVES_r2Dspec_irls(&(spec[i]),&par))
	errormsg("Unknown error returned from UVES_r2Dspec_irls()\n\
\twhen attempting to read in spectrum %d,\n\t%s",i+1,spec[i].file);
      break;
    case FTHIRX:
      if (!UVES_r2Dspec_hirx(&(spec[i]),&par))
	errormsg("Unknown error returned from UVES_r2Dspec_hirx()\n\
\twhen attempting to read in spectrum %d,\n\t%s",i+1,spec[i].file);
      break;
    case FTESOM:
      if (!UVES_r2Dspec_ESOmer(&(spec[i]),&par))
	errormsg("Error returned from UVES_r2Dspec_ESOmer()\n\
\twhen attempting to read in spectrum %d,\n\t%s",i+1,spec[i].file);
      break;
    case FTMAGE:
      if (!UVES_r2Dspec_mage(&(spec[i]),&par))
	errormsg("Error returned from UVES_r2Dspec_mage()\n\
\twhen attempting to read in spectrum %d,\n\t%s",i+1,spec[i].file);
      break;
    case FTIESI:
      if (!UVES_r2Dspec_iresi(&(spec[i]),&par))
	errormsg("Error returned from UVES_r2Dspec_iresi()\n\
\twhen attempting to read in spectrum %d,\n\t%s",i+1,spec[i].file);
      break;
    }
  }
  /* Use object name in first spectrum as object name for combined spectrum */
  if (par.filetype!=FTCOMB) sprintf(cspec.obj,"%s",spec[0].obj);
  else {
    sprintf(cspec.obj,"%s",cspec.abfile);
    if ((cptr=strrchr(cspec.obj,'.'))==NULL) strcat(cspec.obj,"_up");
    else sprintf(cptr,"%s","_up");
  }
  /* If only a combined spectrum is to be read in, read it now and
     ensure certain parameters are set */
  if (par.filetype==FTCOMB) {
    par.combmeth=0;
    if (!UVES_r1Dspec(&cspec,&par))
      errormsg("Error returned from UVES_1Dspec()\n\
\twhen attempting to read in spectrum %d,\n\t%s",i+1,spec[i].file);
  }

  /* If not already known, determine output file names based on object
     name or user-specified prefix */
  if (!strncmp(par.prefix,"\0",1)) {
    if (!strncmp(cspec.UPLfile,"\0",1)) sprintf(cspec.UPLfile,"%s.upl",cspec.obj);
    sprintf(cspec.FITSfile,"%s.fits",cspec.obj);
    sprintf(cspec.DATfile,"%s.dat",cspec.obj);
  } else {
    sprintf(cspec.UPLfile,"%s.upl",par.prefix);
    sprintf(cspec.FITSfile,"%s.fits",par.prefix);
    sprintf(cspec.DATfile,"%s.dat",par.prefix);
  }

  /* If combining ThAr spectra, form an error array from the minimum
     RMS within a sliding window for each order */
  if (par.thar) {
    for (i=0; i<nspec; i++) {
      if (!UVES_thar_sigarray(&(spec[i]),&par))
	errormsg("Error returned from UVES_thar_sigarray()\n\
\twhen operating on spectrum %d, i.e. file\n\t%s",i+1,spec[i].file);
    }
  }

  /* Calculate heliocentric correction for each spectrum */
  if (par.thar<=1) {
    for (i=0; i<nspec; i++) {
      if (!UVES_vhelio(&(spec[i])))
	errormsg("Unknown error returned from UVES_vhelio()");
      // fprintf(stdout,"%s %10.6lf\n",spec[i].abfile,spec[i].vhel);
    }
  }

  /* Median filter the raw error arrays and reject pixels near order
     edges which have grossly underestimated errors */
  if (par.thar<=1) {
    for (i=0; i<nspec; i++) {
      for (j=0; j<spec[i].nor; j++) {
	if (!UVES_order_rejsigedge(&(spec[i].or[j]),&par))
	  errormsg("Error returned from UVES_order_rejsigedge() when\n\
\trejecting bad sigmas from edges of order %d of spectrum %d,\n\ti.e. file\n\t%s",j+1,
		   i+1,spec[i].file);
      }
    }
  }

  /* Sigma-clip each order to reject positive spikes (i.e. cosmic rays) */
  /* Assumption for ThAr spectra is that there are no cosmic rays */
  if (par.thar<=1) {
    for (i=0; i<nspec; i++) {
      for (j=0; j<spec[i].nor; j++) {
	if (!UVES_order_sigclip(&(spec[i].or[j]),&par))
	  errormsg("Error returned from UVES_order_sigclip() when\n\
\tsigma-clipping order %d of spectrum %d,\n\ti.e. file\n\t%s",j+1,i+1,spec[i].file);
      }
    }
  }

  /* Determine wavelength scale of combined spectrum */
  if (nspec) {
    if (!UVES_set_wavelen_scale(spec,nspec,&cspec,&par))
      errormsg("Unknown error while setting wavelength scale of\n\
\tcombined spectrum");
  }

  /* Mask atmospheric/telluric features in observer's reference frame
     if required */
  if (par.atmask) {
    for (i=0; i<nspec; i++) {
      if (!UVES_atmask(&(spec[i]),&amsk,&par))
	errormsg("Error returned from UVES_atmask() when attempting\n\
\tto mask atmospheric features in file\n\t%s",spec[i].file);
    }
  }

  /* Write out original data arrays for all orders to a single FITS file */
  if (par.raw) {
    for (i=0; i<nspec; i++) {
      if (!UVES_wrawFITS(&(spec[i]),&par))
	errormsg("Error returned from UVES_wrawFITS() when attempting\n\
\tto write out original data arrays for spectrum in file\n\t%s",spec[i].file);
    }
  }

  /* Replace ThAr information with synthetic emission line spectrum in
     each un-redispersed order */
  /*
  if (par.thar) {
    for (i=0; i<nspec; i++) {
      if (!UVES_synthThAr(&(spec[i]),&par))
	errormsg("Error returned from UVES_synthThAr() when operating on file\n\
\t%s",spec[i].file);
    }
  }
  */

  /* TEMPORARY FITS TO PIXEL-SPACE ThAr FEATURES */
  /*
  if (par.thar) {
    for (i=0; i<nspec; i++) {
      if (!UVES_tmpfit(&(spec[i]),&par))
	errormsg("Error returned from UVES_tmpfit() when operating on file\n\t%s",
		 spec[i].file);
    }
  }
  */

  /* Redispers each order according to combined spectrum's wavelength scale */
  for (i=0; i<nspec; i++) {
    if (!UVES_redispers(&(spec[i]),&cspec,spec[0].distort_seed,&par))
      errormsg("Unknown error returned from UVES_redispers()\n\
\twhen operating on file\n\t%s",spec[i].file);
  }

  /* TEMPORARY: Write out start and end wavelengths of each order */
  /*
  for (i=0; i<nspec; i++) {
    for (j=0; j<spec[i].nor; j++) {
      fprintf(stdout,"%s  order=%2d   starl_wl=%.4lf   end_wl=%.4lf\n",
	      spec[i].abfile,j+1,cspec.wl[spec[i].or[j].csidx],
	      cspec.wl[spec[i].or[j].ceidx]);
    }
  }
  */

  /* Do some simple statistics on each order, including calculating a
     median error array which will be used later for combining the
     different spectra */
  for (i=0; i<nspec; i++) {
    for (j=0; j<spec[i].nor; j++) {
      if (!UVES_order_stats(&(spec[i].or[j]),&par))
	errormsg("Failed to calculate statistics for order\n\
\t%d in spectrum %d, i.e. file\n\t%s",j+1,i+1,spec[i].file);
    }
  }

  /* Merge the orders (in a simple step-function way) of the ThAr exposures */
  if (par.thar==1) {
    for (i=0; i<nspec; i++) {
      if (!UVES_merge_thar(&(spec[i]),&cspec,&par))
	errormsg("Error returned from UVES_merge_thar when operating\n\
\ton spectrum %d in file\n\t%s",i+1,spec[i].file);
    }
  }

  /* Combine spectra in either of two different ways */
  if (!par.combmeth) {
    /* Combine the orders and exposures without first finding a continuum */
    if (!UVES_combine_region(spec,nspec,&cspec,&par))
      errormsg("Error returned from UVES_combine_region()");
    /* Set a continuum for the combined spectrum */
    if (!UVES_cspec_cont(spec,nspec,&cspec,&par))
      errormsg("Error returned from UVES_cspec_cont()");
  } else {
    /* Make a rough continuum for each order of every spectrum */
    if (!UVES_order_cont(spec,nspec,&par))
      errormsg("Unknown error returned from UVES_order_cont()");
    /* Make a unified continuum for the combined spectrum */
    if (!UVES_combine_cont(spec,nspec,&cspec,0,&par))
      errormsg("Unknown error returned from UVES_combine_cont()");
    /* Combined all orders and spectra together */
    if (!UVES_combine_spec(spec,nspec,&cspec,&par))
      errormsg("Problem combining spectra in UVES_combine_spec()");
  }

  /* Carry out historical actions from previous reduction attempts */
  if (nact && !par.replay) {
    if (!UVES_past_actions(spec,nspec,&cspec,act,nact,0,&par))
      errormsg("Unknown error returned from UVES_past_actions()");
  }
  nact_save=nact;

  /* Replace ThAr information with synthetic emission line spectrum */
  /*
  if (par.thar) {
    for (i=0; i<nspec; i++) {
      if (!UVES_combsynthThAr(&cspec,&par))
	errormsg("Error returned from UVES_combsynthThAr()\n\
\twhen operating on file\n\t%s",spec[i].file);
    }
  }
  */

  /* Skip display if user only wants to save UPL and FITS output */
  if (par.save!=2) {
    /* Initialize plotting */
    UVES_pgenv_init(&plenv,&cp);
    sprintf(cp.ptitle,"%s: Spectrum Navigator",progname);
    pg_open(&plenv,"?\0",cp.ptitle,0); cpgqid(&(cp.pgid)); rp.init=rp.cact=rp.exit=0;

    /* Plot spectrum and allow user to perform certain actions on spectrum */
    while (!cp.exit) {
      if (nact && par.replay) {
        if (!UVES_replay_control(spec,nspec,&cspec,&cp,&rp,&act,&nact,&nact_save,
				 &par))
	  errormsg("Problem controlling interactive action replay");
	if (rp.exit==1) par.replay=0;
      }
      if (!UVES_plot_cspec(spec,nspec,&cspec,&cp,&rp,&act,&nact,&nact_save,&par))
	errormsg("Problem plotting combined spectrum");
      if (cp.refre) cp.refre=cp.exit=0;
      else if (cp.rscsp && !par.combmeth) {
	/* Auto-rescale orders after 1 or more have been manually rescaled */
	/* Determine which action were are currently dealing with */
	if (par.replay && rp.exit==2) cact=rp.cact-1;
	else cact=nact-1;
	/* Rescale the orders */
	status=UVES_rescale_region(spec,nspec,&cspec,act,cact+1,cp.scalclip,
				   cp.scalerr,&par);
	switch (status) {
	case 0: errormsg("Error returned from UVES_rescale_region()"); break;
	case 2: warnmsg("No manually rescaled orders since\n\
\tlast recombination"); break;
	case 1:
	  /* Remember autorescaling action by filling in an action report */
	  if (!(act=(action *)realloc(act,(size_t)((nact+1)*sizeof(action)))))
	    errormsg("Cannot increase memory for act array\n\
\tto size %d",nact+1);
	  /* If in paused replay mode then shift later actions to end */
	  if (par.replay && rp.exit==2) {
	    for (j=nact-1; j>=rp.cact; j--) act[j+1]=act[j];
	    cact=rp.cact++; nact++;
	  } else cact=nact++;
	  act[cact].act=ARACT; act[cact].nxyp=0;
	  act[cact].d[0]=cp.scalclip; act[cact].d[1]=cp.scalerr;
	  act[cact].rcmb=act[cact].val=act[cact].nordact=1; break;
	}
	cp.rscsp=cp.exit=0;
      } else if (cp.rccsp) {
	/* Recombine all orders and spectra */
	if (!UVES_combine_spec(spec,nspec,&cspec,&par))
	  errormsg("Problem combining spectra in UVES_combine_spec()");
	cp.rccsp=cp.rscsp=cp.exit=0;
      }
    }
  
    /* Close plotting */
    cpgclos();
  }

  /* Write out log file to record file list and all historical actions made */
  if (par.save) {
    if (par.save==2) par.save=0;
    if (!UVES_wUPLfile(cspec.UPLfile,spec,nspec,act,nact,&cspec,&par))
      errormsg("Unknown error returned from UVES_wUPLfile()");
    if (!UVES_wFITSfile(spec,nspec,&cspec,&par))
      errormsg("Unknown error returned from UVES_wFITSfile()");
    if (par.dat && !UVES_wDATfile(cspec.DATfile,&cspec))
      errormsg("Unknown error returned from UVES_wDATfile()");
  } else if (nact!=nact_save) {
    sprintf(query,"y"); get_input("Save UPL and 1D FITS file?","%s",query);
    if (strncmp(query,"n",1)) {
      if (!UVES_wUPLfile(cspec.UPLfile,spec,nspec,act,nact,&cspec,&par))
	errormsg("Unknown error returned from UVES_wUPLfile()");
      if (!UVES_wFITSfile(spec,nspec,&cspec,&par))
	errormsg("Unknown error returned from UVES_wFITSfile()");
      if (par.dat && !UVES_wDATfile(cspec.DATfile,&cspec))
	errormsg("Unknown error returned from UVES_wDATfile()");
    }
  }

  return 1;

}
