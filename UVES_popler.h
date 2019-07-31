/***************************************************************************
* Definitions, structures and function prototypes for UVES_popler
***************************************************************************/

/* INCLUDE FILES */
#include "pg_plot.h"

/* DEFINITIONS */
/* Functions */
#define NINT(a)       (int)(0.5+a)
#define MAX(a,b)      a>b ? a : b
#define MIN(a,b)      a<b ? a : b
#define RMDR(a,b)     a-b*(a/b)
#define PMOD(a,b,c,d) (a+b)>c ? d : (a+b)
#define DSWAP(a,b)    dtmp=(a);(a)=(b);(b)=dtmp
#define FSWAP(a,b)    ftmp=(a);(a)=(b);(b)=ftmp
#define ISWAP(a,b)    itmp=(a);(a)=(b);(b)=itmp

/* Version number, date created and UVES_popler website */
#define VERSION       1.04    /* Version number */
#define DATECREATE    "16 June 2019"
#define WWW           "www.astronomy.swin.edu.au/~mmurphy/UVES_popler.html"

/* General (non-option) Parameters */
#define DRNDTOL       1.e-15  /* Rounding tolerance for double precision     */
#define INFIN         1.e32   /* Effective infinity parameter                */
#define MINUSE      100       /* Minimum useful # of pixels in an order      */
#define UTTOL         4.0     /* Tol. for UT diffs in header vals [min]      */
#define WLTOL         1.e-10  /* Tol. for wavelength scale linearity         */
#define LYA        1215.67    /* Lyman-alpha wavelength [A]                  */
#define NWPOL        12       /* Number of wavelength polynomial coefficients*/
#define NINTP         4       /* Num.pnts.poly.interp. 4 invert.wavelen.poly.*/
                              /*    MUST BE >1 AND <=NWPOL                   */
#define TNCOLCOM       7      /* Num. cols in comb. spec. info. output table */ 
#define TNCOLWCO       2      /* Num. cols in wav. coverage info. output tab.*/ 
#define TNCOLEXP      55      /* Num. cols in exposure info. output table    */ 
#define UVESBORR    500.0     /* Cross-over wave. (nm) for blue/red settings */
#define UVESPIXSCALB  0.215   /* Central order pixel scale of blue UVES arm  */
#define UVESPIXSCALR  0.155   /* Central order pixel scale of red UVES arm   */
#define HARPSNSPATPIX 4.0     /* Eff. # spatial pix per HARPS-S extracted    */
                              /*    spectral pix                             */
#define HARPNNSPATPIX 4.0     /* Eff. # spatial pix per HARPS-N extracted    */
                              /*    spectral pix                             */
#define WAVGRIDSIZE 1000.0    /* Size of coarse wavelength grid for measuring*/
                              /*    statistics on combined spectrum [km/s]   */
#define WAVGRIDSTART 3000.0   /* Starting wavelength of coarse grid [A]      */
#define WAVGRIDEND  11000.0   /* Ending wavelength of coarse grid [A]        */
#define WAVCOVMINREG 100.0    /* Min. size of chunk of wav. cov. map         */
#define WAVCOVMAXGAP  10.0    /* Max. gap size within chunk in wav. cov. map */

/* Command line options */
 /* General spectrum specifications */
#define COMBMETH      0       /* Switch to combine exposures using order-fits*/
#define LINEAR        0       /* Log-lin (0) or linear (1) dispersion switch */
#define FILETYPE      0       /* File origin: UVES(0), IRAF(1), MAKEE(2),    */
                              /* IRAFLSS(3), HIREDUX(4), ESOMERGED(5),MAGE(6)*/
                              /* KODIAQ(7), IRAFESI (8), COMB (10), mixed(-1)*/
#define HELIO         0       /* Files in obs. (0) or heliocen. (1) frame    */
#define VACWL         0       /* Files in air (0) or vac. (1) wavelengths    */
#define THAR          0       /* Swith to combine different ThAr exposures   */
#define ZEM           0.0     /* Emission redshift of QSO                    */
 /* Initial sigma-clipping of orders (not done when using -thar) */
#define NORDCLIP      3       /* Most positive pixels removed before rms     */
                              /*    determined                               */
#define NORDSIG      35       /* Window width for finding rms -              */
                              /*    MUST BE ODD AND MUST BE > 5              */
#define ORDSIG        5.0     /* Clip level, units of the rms in window      */
#define ORDSIGNBR     3.0     /* Clip level for neighbouring points          */
#define ORDSIGZERO    3.0     /* Only clip when this many sigma above zero   */
#define ORDMEDFRAC    0.10    /* Median filter width relative to order length*/
#define ORDMEDREJ     0.05    /* Reject pixels this much smaller than median */
 /* Spectrum/order combination parameters */
#define CLIPSIG       3.0     /* Clip level for sigma clip                   */
#define LRGERR        4.0     /* Rejection ratio for large errors            */
#define RANKSPEC      0       /* For combmeth=0, 0=Normal scaling algorithm; */
                              /*    1=don't scale high SNR ord of first spec */
#define SCALMETH     1        /* Switch to use full chisq. min. order scaling*/
#define SCALCLIP     3.00     /* Clip level for inter-order scale-factor     */
#define SCALERR      0.20     /* Lim. for rel. err. on inter-order scale-fac.*/
#define NSCALCLIP    5        /* Max. rej. scaling iterations (scalmeth=1)   */
 /* Order continuum fitting method parameters */
#define CONTPCTL      0.10    /* Lower %ile removed for init. ord. cont. fit.*/
#define CONTSIGL      1.4     /* Lower rej. level for order cont. fitting    */
#define CONTSIGU      3.0     /* Upper rej. level for order cont. fitting    */
#define CONTFTYP      2       /* Type of fit (1=Poly, 2=Cheb, 3=Lege)        */
#define CONTORD       4       /* Order of fitted continuum (1=const.)        */
#define CONTSMOOTH  501       /* No. smoothing pixels for joining order cont.*/
#define CONTWGT     500       /* No. pixels in weighting ends of orders      */
 /* Combined spectrum continuum fitting method parameters */
#define CORDLYA       4       /* Continuum fitting order for Lya cont.       */
#define CORDRED       6       /* Continuum fitting order for red cont.       */
#define FTYPLYA       2       /* Type of fit below Lya (1=Pol, 2=Che, 3=Leg) */
#define FTYPRED       2       /* Type of fit above Lya (1=Pol, 2=Che, 3=Leg) */
#define NOCONT        0       /* 1=Do not attempt auto-cont. fit, 0=default  */
#define PCTLLYA       0.30    /* Lower %ile removed for initial Lya cont. fit*/
#define PCTLRED       0.10    /* Lower %ile removed for initial red cont. fit*/
#define RSIGLYAL      1.2     /* Lower rej. limit [sigma] for Lya cont.      */
#define RSIGLYAU      3.0     /* Upper rej. limit [sigma] for Lya cont.      */
#define RSIGREDL      1.4     /* Lower rej. limit [sigma] for red cont.      */
#define RSIGREDU      3.0     /* Upper rej. limit [sigma] for red cont.      */
#define VCLYA     20000.0     /* Vel. scale for continuum below Lya [kms/s]  */
#define VCRED      2500.0     /* Vel. scale for continuum above Lya [kms/s]  */
#define VLYA       5000.0     /* Vel. below Lya where red cont. begins[kms/s]*/
 /* Data file options */
#define DAT           0       /* Option for writing out ASCII data file      */
#define RAW           0       /* Option for writing out raw order data       */
#define SAVE          0       /* Option for forced save of UPL and FITS file */
/* Mode options */
#define REPLAY        0       /* Option for entering replay mode             */
#define MACMAP        0       /* Use a Macmap file (1) or not (0)            */
#define MACMAPFILE "UVES_headsort.macmap"
                              /* Default name for Macmap from UVES_headsort  */
#define VSHIFT        0       /* Apply velocity shifts from a to spectra     */
#define VSHIFTFILE "UVES_popler.vshift"
                              /* Default name for velocity shift file        */
#define SCALE         0       /* Apply flux/error scaling to spectra         */
#define SCALEFILE "UVES_popler.scale"
                              /* Default name for scalings file              */
#define ATMASK        0       /* Apply mask to atmospheric features          */
#define ATMASKFILE "UVES_popler.atmomask"
                              /* Default name for atmo. mask file            */
#define ATMASKMODE    0       /* Default mode for atmo. masking = after comb.*/
#define DISTORT       0       /* Option for applying, or removing existing   */
                              /* distortions to wavelength scales of orders  */

/* Size of sliding median window for calculating median error array */
#define NMEDERR      13

/* General continuum fitting parameters */
#define MAXITERFIT  100    

/* Interactive continuum fitting parameters */
#define FFRAC         0.10    /* Fraction of plotted order to disregard when */
                              /*   forming first guesses for fit region      */
#define FBFRAC        0.20    /* Fraction of fit region to use aa blending   */
                              /*   region with old continuum                 */
#define SSCALE        0.01    /* Fraction by which to minimally scale order  */ 
#define BSCALE        0.10    /* Fraction by which to maximally scale order  */ 

/* Error array labelling */
#define RCLIP        -1      /* Clipped during read from original file       */
#define OCLIP        -2      /* Clipped during single order sigma-clip       */
#define ACLIP        -3      /* Masked atmospheric/telluric features         */
#define LCLIP        -4      /* Clipped due to comparitively large error     */
#define NCLIP        -5      /* No contributing pixels                       */
#define SCLIP        -6      /* Clipped during combination sigma-clip        */
#define CCLIP        -7      /* Cosmetically clipped                         */
#define PCLIP        -8      /* Clipped when refitting continuum after plot  */

/* Action labelling */
#define COACT         1      /* Clip pixels from order                       */
#define CCACT         2      /* Clip pixels from combined spectrum           */
#define UOACT         3      /* Unclip pixels from order                     */
#define UCACT         4      /* Unclip pixels from combined spectrum         */
#define FOACT         5      /* Fit new continuum to part of order           */
#define FCACT         6      /* Fit new continuum 2 part of combined spectrum*/
#define IOACT         7      /* Interpolate new continuum to part of order   */
#define ICACT         8      /* Interpolate new cont. 2 part of comb. spec.  */
#define SOACT         9      /* Scale order flux                             */
#define ARACT        10      /* Autorescale orders after user has scaled some*/
#define NCACT        11      /* Create new continuum                         */
#define AAACT        12      /* Switch alternative arrays to main arrays     */

/* Filetype labelling */
#define FTMIX        -1      /* A mix of different filetypes is to be used   */
#define FTUVES        0      /* UVES pipeline products                       */
#define FTIRAF        1      /* IRAF reduction products                      */
#define FTMAKE        2      /* MAKEE pipeline products                      */
#define FTIRLS        3      /* IRAF LSS reduction products                  */
#define FTHIRX        4      /* HIREDUX reduction products                   */
#define FTESOM        5      /* ESOMERGE reduction products                  */
#define FTKODI        6      /* KODIAQ products                              */
#define FTMAGE        7      /* MAGE pipeline products                       */
#define FTIESI        8      /* IRAF ESI pipeline products                   */
#define FTHARP        9      /* HARPS (S & N) pipeline products              */
#define FTESPR       10      /* ESPRESSO pipeline products                   */
#define FTCOMB       11      /* Combined spectrum                            */

/* Fit labelling */
#define FITPOL        1      /* Polynomial fit                               */
#define FITCHE        2      /* Chebyshev fit                                */
#define FITLEG        3      /* Legendre fit                                 */

/* Continuum options */
#define CNTNON        1      /* User says don't fit auto-continuum           */
#define CNTSRD        2      /* No auto-cont. fitted, cont. read in from file*/

/* General plotting parameters */
#define W_WIDTH     16.0
#define W_ASP        0.75
#define C_H          1.4
#define L_W          2.0
#define VPU          0.96
#define VPD          0.10
#define VPL          0.07
#define VPR          0.99
#define N_X_SUB      1
#define N_Y_SUB      1

/* Combined plot default parameters */
#define VPD0        -0.02    /* Lower view-port limit for panel 0       */
#define VPU0         0.04    /* Upper view-port limit for panel 0       */
#define VPD1         0.10    /* Lower view-port limit for panel 1       */
#define VPU1         0.22    /* Upper view-port limit for panel 1       */
#define VPD2         0.22    /* Lower view-port limit for panel 2       */
#define VPU2         0.34    /* Upper view-port limit for panel 2       */
#define VPD3         0.34    /* Lower view-port limit for panel 3       */
#define VPU3         0.78    /* Upper view-port limit for panel 3       */
#define VPD4         0.78    /* Lower view-port limit for panel 4       */
#define VPU4         1.00    /* Upper view-port limit for panel 4       */
#define SFAC        10.00    /* Degredation factor for low-res spectrum */
#define CP_NP     5000       /* Number of pixels in detailed spec.      */

/* Sub-spectrum selection panel parameters */
#define NSPPAGE     25       /* Number of spectra listed per page       */
#define SPSIZE       0.75    /* Ratio of widths of selection panel and  */
                             /*    combined plot                        */
#define RPSIZE       0.15    /* Ratio of widths of replay plot and      */
                             /*    combined plot                        */

/* STRUCTURES */
typedef struct Params {
  double    clipsig;              /* Rejection level in sig for spec. comb.n.*/
  double    contpctl;             /* Lower %ile to remove for intial cont fit*/
  double    contsigu;             /* Upper rej. level in sig for cont fit    */
  double    contsigl;             /* Lower rej. level in sig for cont fit    */
  double    disp;                 /* Dispersion in either km/s or Angstroms  */
  double    lrgerr;               /* Rejtn ratio for large errs during comb. */
  double    ordsig;               /* Clip level for orders - units of rms    */
  double    ordsignbr;            /* Clip level for neighbouring pixels      */
  double    ordsigzero;           /* Only clip when this many sig above zero */
  double    ordmedfrac;           /* Med, filt. width rel. to order length   */
  double    ordmedrej;            /* Reject pixels this much smaller than med*/
  double    pctllya;              /* Low. %ile to remove for initial Lya fit */
  double    pctlred;              /* Low. %ile to remove for initial red fit */
  double    rsiglyau;             /* Upper rej. limit [sigma] for Lya cont.  */
  double    rsiglyal;             /* Lower rej. limit [sigma] for Lya cont.  */
  double    rsigredu;             /* Upper rej. limit [sigma] for red cont.  */
  double    rsigredl;             /* Lower rej. limit [sigma] for red cont.  */
  double    scalclip;             /* Clip level on inter-order scale-factor  */
  double    scalerr;              /* Rel. err. lim. on inter-order scale-fac.*/
  double    vclya;                /* Velocity scale for continuum below Lya  */
  double    vcred;                /* Velocity scale for continuum below red  */
  double    vlya;                 /* Vel. bel. Lya where red cont. pars apply*/
  double    zem;                  /* Emission redshift of QSO                */
  double    version;              /* Version of UPL file read in             */
  int       atmask;               /* Use an atmospheric features input file  */
  int       combmeth;             /* Switch to use order-fitting method      */
  int       contftyp;             /* Type of cont. fit (1=pol, 2=che, 3=leg) */
  int       contord;              /* Order of fitted continuum               */
  int       contwgt;              /* Weighting scale for ends of orders      */
  int       cordlya;              /* Continuum fitting order for Lya cont.   */
  int       cordred;              /* Continuum fitting order for red cont.   */
  int       dat;                  /* Toggle for writing of ASCII data file   */
  int       distort;              /* Toggle for random wavelength distortions*/
  int       filetype;             /* File origin: UVES (0) IRAF (1) MAKEE (2)*/
                                  /*   COMB (9) mixed (-1)                   */
  int       ftyplya;              /* Fit type below Lya (1=pol, 2=che, 3=leg)*/
  int       ftypred;              /* Fit type above Lya (1=pol, 2=che, 3=leg)*/
  int       helio;                /* Files in obs. (0) or heliocen. (1) frame*/
  int       linear;               /* Linear (1) or log-lin (0) dispersion    */
  int       macmap;               /* Use a Macmap input file (1) or not (0)  */
  int       vshift;               /* Use a velocity shift input file? 1 or 0 */
  int       nocont;               /* Cont. fit (0) or no auto-cont. fit (1)  */
  int       nordclip;             /* No. +ve pixels remvd before rms found   */
  int       nordsig;              /* Window width for sigma-clipping orders  */
  int       nscalclip;            /* Max. rej. iters for scaling (scalmeth=1)*/
  int       raw;                  /* Option for writing out raw order data   */
  int       rankspec;             /* Switch between scaling order choice     */
  int       save;                 /* Force UPL and FITS to save upon quitting*/
  int       scale;                /* Switch to read in flux/error scale file */
  int       scalmeth;             /* Switch to use full chisq. min. scaling  */
  int       thar;                 /* Allow combination of ThAr cal. exps.    */
  int       replay;               /* Replay mode?                            */
  int       vacwl;                /* Files in air (0) or vacuum (1) waves    */
  int       backvers;             /* Flag (1) for badging output UPL with    */
				  /*   input UPL file's version number       */
  char      atmaskfile[NAMELEN];  /* Name of atmospheric feature mask file   */
  char      prefix[NAMELEN];      /* Prefix for output file names            */
  char      macmapfile[NAMELEN];  /* Name of Macmap file                     */
  char      scalefile[NAMELEN];   /* Name of flux/error scaling file         */
  char      vshiftfile[NAMELEN];  /* Name of velocity shift file             */
} params;

typedef struct EchOrder {
  double   fvhlwl;                /* Left edge of first pixel in vac-helio   */
  double   medsnr;                /* Median S/N of valid, redispersed pixels */
  double   inscl;                 /* Input scaling factor from input file    */
  double   scl;                   /* Scaling factor applied to order         */
  double   oscl;                  /* Old scaling factor before last action   */
  double   seeing;                /* FWHM of spatial prof. in object extractn*/
  double   wpol[NWPOL];           /* Wavelength polynomial coefficients      */
  double   *wl;                   /* Pntr to raw wavelength array read from  */
                                  /*  input files for HIRES REDUX, KODIAQ,   */
                                  /*  MAGE, ESPRESSO only                    */
  double   *vhwl;                 /* Pntr to vac-hel arr of raw spec         */
  double   *vhrwl;                /* Pntr to right edge vac-hel array        */
  double   *fl;                   /* Pointer to flux array of raw spec       */
  double   *rdfl;                 /* Pointer to redispersed flux array       */
  double   *er;                   /* Pointer to error array of raw spec      */
  double   *rder;                 /* Pointer to redispersed error array      */
  double   *res;                  /* Pointer to resolution array             */
  double   *rdres;                /* Pointer to redispersed resolution array */
  double   *rdef;                 /* Pointer to redisp. expected fluc. array */
  double   *rdco;                 /* Pointer to redispersed continuum array  */
  double   *ordco;                /* Old redispersed cont. array b4 last act.*/
  double   *rdme;                 /* Pointer to redispersed median error arr */
  double   *th;                   /* Pointer to ThAr flux array of raw spec  */
  double   *ter;                  /* Pointer to ThAr error array of raw spec */
  double   *rdth;                 /* Pointer to redispersed ThAr flux array  */
  double   *rdter;                /* Pointer to redispersed ThAr error array */
  double   *rdfwhm;               /* Pointer to redisp. resol. (FWHM) array  */
                                  /** Arrays for pre-v0.74 backwards compat **/
  double   *rdfl_a;               /* Alternative redispersed flux array      */
  double   *rder_a;               /* Alternative redispersed error array     */
  double   *rdef_a;               /* Alternative redisp. expected fluc. array*/
  double   *rdme_a;               /* Alternative redispersed median error arr*/
  int      id;                    /* Identity number for this order          */
  int      sid;                   /* ID of spec. to which this order belongs */
  int      insclbefaft;           /* Input scaling: 0=none, 1=before, 2=after*/
                                  /*    automatic and manual actions         */
  int      insclfe;               /* Input scaling: 0=flux, 1=error, 2=both  */
  int      ascl;                  /* Flag whether order was auto-rescaled    */
  int      idxoff;                /* Offset in pixel index when reading FITS */
                                  /*    files (=CRVAL1 value, usuaully)      */
  int      np;                    /* Number of pixels in raw arrays          */
  int      nrdp;                  /* Number of pixels in redispersed arrays  */
  int      nuse;                  /* Useful (i.e. er!=0.0) pixels in order   */
  int      nwpol;                 /* Number of wavelength polynomial coeffs  */
  int      crank;                 /* Combination rank of order for combmeth=0*/
  int      csidx;                 /* Index of first pixel in cspec to which  */
                                  /*    pixels in this order contribute      */
  int      ceidx;                 /* Index of last pixel in cspec to which   */
                                  /*    pixels in this order contribute      */
  int      revdisp;               /* Dispersion positive(0) or negative(1)?  */
  int      *st;                   /* Status of each pixel                    */
  int      *rdst;                 /* Status of each redispersed pixel        */
  int      *ordst;                /* Old status before last action           */
  int      *tst;                  /* Status of each ThAr pixel               */
  int      *rdtst;                /* Status of each redispersed ThAr pixel   */
} echorder;

typedef struct SkySpec {
  int      nor;                   /* Number of orders                        */
  char     path[LNGSTRLEN];       /* Full path of sky flux file              */
  char     file[LNGSTRLEN];       /* Name of flux file                       */
  char     abfile[LNGSTRLEN];     /* Name of flux file without full path     */
  echorder *or;                   /* Pointer to echelle order array          */
} skyspec;

typedef struct ThArSpec {
  double   fwl;                   /* Vac. wavelength of first pixel's center */
  double   lwl;                   /* Vac. wavelength of last pixel's center  */
  double   *fl;                   /* Pointer to flux array                   */
  double   *er;                   /* Pointer to error array                  */
  int      np;                    /* Number of pixels in arrays              */
  int      csidx;                 /* Index of first pixel in cspec for first */
                                  /*    ThAr spec. pixel                     */
  int      ceidx;                 /* Index of last pixel in cspec for last   */
                                  /*    ThAr spec. pixel                     */
  int      *st;                   /* Status of each pixel                    */
} tharspec;

typedef struct ThArSet {
  double   *x;                    /* Fitted pixel position                   */
  double   *w;                    /* Fitted FWHM width of line [pix]         */
  double   *dis;                  /* Dispersion [A/pix]                      */
  double   *wlf;                  /* Wavelength from fit [A]                 */
  double   *resid;                /* Residual between fit and lab wavelen [A]*/
  double   *wlsf;                 /* FWHM of line-spread function            */
  int      binx;                  /* CCD Binning in spectral direction       */
  int      n;                     /* Number of lines in set                  */
  int      np;                    /* Number of lines in polynomial solution  */
  int      *stp;                  /* Status as line used in polynomial fit   */
} tharset;

typedef struct Spectrum {
  double   jd;                    /* Julian day, start of exp. [days]        */
  double   ut;                    /* UT from header                          */
  double   cwl;                   /* Nominal central wavelength of setting   */
  double   sw;                    /* Slit-width in arcsec                    */
  double   temp;                  /* Temperature at echelle grating          */
  double   pres;                  /* Pressure inside spectrograph            */
  double   pixscal;               /* Pixel scale in spatial direction ["/pix]*/
  double   etime;                 /* Exposure time [s]                       */
  double   ra;                    /* Right ascension of object [hrs]         */
  double   dec;                   /* Declination of object [deg]             */
  double   epoch;                 /* Epoch of observation start              */
  double   equ;                   /* Equinox of RA and DEC                   */
  double   lat;                   /* Latitude of Observatory [+=North, deg]  */
  double   lon;                   /* Longitude of Observatory [+=West, deg]  */
  double   alt;                   /* Altitude of Observatory (meters)        */
  double   vhel;                  /* Heliocen. vel. at middle of exp. (calc) */
  double   vhel_head;             /* Heliocen. vel. from header              */
  double   vbar_head;             /* Barycen. vel. from header               */
  double   vshift;                /* Velocity shift applied by user [km/s].  */
  double   vslope;                /* Velocity slope " [km/s/1000Ang].        */
  double   refwav;                /* Ref. wav. for vel. slope [Ang].         */
  double   airmass[2];            /* Airmass at start [0] and end [1] of exp.*/  
  double   seeing[3];             /* Seeing [arcsec]: [0]=start(header),[1]= */
                                  /* end(header),[3]=med. extracted          */  
  double   moon_ra;               /* RA of Moon [hrs]                        */
  double   moon_dec;              /* Dec. of Moon [deg]                      */
  double   moon_ang;              /* Angle between object and Moon           */
  double   moon_phase;            /* Moon phase as fraction of period        */
  double   arcfwhm;               /* Median FWHM of ThAr arc lines [km/s]    */
  double   wc_resid;              /* Median residual wavelen. solution [m/s] */
  double   wc_jd;                 /* Wavcal Julian day, start of exp. [days] */
  double   wc_ut;                 /* Wavcal UT start of exposure [days]      */
  double   wc_etime;              /* Wavcal exposure time [s]                */
  double   wc_sw;                 /* Wavcal slit-width in arcsec             */
  double   wc_temp;               /* Temperature at echelle grating for ThAr */
  double   wc_pres;               /* Pressure inside spectrograph for ThAr   */
  int      id;                    /* Identity number for this spectrum       */
  int      comb;                  /* Is this spectrum to be combined?        */
  int      ftype;                 /* File type: UVES=0, IRAF=1, MAKEE=2,     */
                                  /*   COMB=9, mixed=-1                      */
  int      fvers;                 /* File version (e.g. for UVES, this is the*/
                                  /*   pipeline version which reduced it)    */
  int      binx;                  /* CCD Binning in spectral direction       */
  int      biny;                  /* CCD Binning in spatial direction        */
  int      year;                  /* Year of observation start (from mjd)    */
  int      month;                 /* Month of observation start (from mjd)   */
  int      day;                   /* Day of observation start (from mjd)     */
  int      nor;                   /* Number of orders                        */
  int      nhead_ori;             /* Num. keys in main header of flux file   */
  int      skysub;                /* Flag to indicate need for separate sky- */
                                  /*   subtraction (=1).                     */
  int      obid;                  /* Observation block ID (for UVES files)   */
  int      encoder;               /* Grating encoder value                   */
  int      wc_year;               /* Wavcal year of obs. start (from mjd)    */
  int      wc_month;              /* Wavcal month of obs. start (from mjd)   */
  int      wc_day;                /* Wavcal day of obs. start (from mjd)     */
  int      wc_encoder;            /* Grating encoder value for wavcal exp.   */
  long     distort_seed;          /* Random number seed for wavelength       */
                                  /*   distortions                           */
  char     path[LNGSTRLEN];       /* Full path of flux and error file        */
  char     file[LNGSTRLEN];       /* Name of flux file                       */
  char     abfile[LNGSTRLEN];     /* Name of flux file without full path     */
  char     erfile[LNGSTRLEN];     /* Name of error array file                */
  char     aberfile[LNGSTRLEN];   /* Name of error file without full path    */
  char     thfile[LNGSTRLEN];     /* Name of ThAr file                       */
  char     abthfile[LNGSTRLEN];   /* Name of ThAr file without full path     */
  char     wlfile[LNGSTRLEN];     /* Name of wavelength solution array file  */
  char     abwlfile[LNGSTRLEN];   /* Name of wavelen. file without full path */
  char     arfile[LNGSTRLEN];     /* Name of archived file                   */
  char     wlarfile[LNGSTRLEN];   /* Name of wavelength sol. archived file   */
  char     tharfile[LNGSTRLEN];   /* Name of ThAr archived file              */
  char     obj[NAMELEN];          /* Object name                             */
  char     progid[NAMELEN];       /* Program ID                              */
  char     dprtech[NAMELEN];      /* DPR.TECH                                */
  char     dprtype[NAMELEN];      /* DPR.TYPE                                */
  char     dprcatg[NAMELEN];      /* DPR.CATG                                */
  char     inspath[NAMELEN];      /* Instrument light path used              */
  char     insmode[NAMELEN];      /* Instrument mode used (e.g. dichr#2)     */
  char     insdrot[NAMELEN];      /* Instrument derotator mode used          */
  char     inscl[VHUGESTRLEN];    /* Flux/error scaling string for all orders*/
  char     **head_ori;            /* Original flux file main header          */
  skyspec  ss;                    /* Sky spectrum                            */
  tharspec th;                    /* Order-merged 1D ThAr spectrum           */
  tharset  ts;                    /* Set of ThAr line parameters             */
  echorder *or;                   /* Pointer to echelle order array          */
} spectrum;

typedef struct CSpectrum {
  double   dv;                    /* Dispersion in velocity space            */
  double   dwl;                   /* Dispersion in wavelength space          */
  double   flwl;                  /* First left edge vac-helio wavelen       */
  double   texp;                  /* Total exposure time                     */
  double   *wl;                   /* Pointer to vac-helio wavelength array   */
  double   *rwl;                  /* Pointer to right edge vac-helio array   */
  double   *fl;                   /* Pointer to flux array                   */
  double   *er;                   /* Pointer to error array                  */
  double   *ef;                   /* Pointer to expected pixel-to-pixel      */
                                  /*    fluctuation array                    */
  double   *co;                   /* Pointer to continuum array              */
  double   *res;                  /* Pointer to resolution array             */
  double   *oco;                  /* Old continuum array before last action  */
  double   *csq;                  /* Chisq. array bef. sig-clip during comb. */
  double   *ccsq;                 /* Chisq. array aft. sig-clip during comb. */
  double   *no;                   /* Pointer to normalized flux array        */
  double   *ne;                   /* Pointer to normalized error array       */
  double   *nf;                   /* Pointer to normalized fluctuation array */
  double   *cwav;                 /* Coarse wavelength array (for stats)     */
  double   *csnrmax;              /* Coarse maxium SNR array (for stats)     */
  double   *csnrmed;              /* Coarse median SNR array (for stats)     */
  double   *ccnrmax;              /* Coarse maxium CNR array (for stats)     */
  double   *ccnrmed;              /* Coarse median CNR array (for stats)     */
  double   *cresnom;              /* Coarse nominal resolution array (for sta*/
  double   *cresarc;              /* Coarse arc resolution array (for stats) */
  double   **wavcov;              /* Wavelength coverage map                 */
  int      np;                    /* Number of pixels in arrays              */
  int      nc;                    /* Number of chunks in coarse stats arrays */
  int      nwavcov;               /* Number of regions in wav. cov. matrix   */
  int      nexp;                  /* Number of exposures (not files)         */
  int      *ncb;                  /* Array of no. of contribut. pxls to comb.*/
  int      *nccb;                 /* No. contribut. pxls to comb. after clip */
  int      *st;                   /* Status of each pixel                    */
  int      *ost;                  /* Old status before last action           */
  char     path[LNGSTRLEN];       /* Path for file for filetype=FTCOMB       */
  char     file[LNGSTRLEN];       /* Original filename for filetype=FTCOMB   */
  char     abfile[LNGSTRLEN];     /* Same as above but without full path     */
  char     obj[NAMELEN];          /* Object name                             */
  char     DATfile[LNGSTRLEN];    /* UVES_popler DAT file name               */
  char     FITSfile[LNGSTRLEN];   /* UVES_popler FITS file name              */
  char     UPLfile[LNGSTRLEN];    /* UVES_popler log file name               */
} cspectrum;

typedef struct TwoXYP {
  double x1;                      /* x-position 1                            */
  double y1;                      /* y-position 1                            */
  double x2;                      /* x-position 2                            */
  double y2;                      /* y-position 2                            */
  int    i;                       /* Integer describing what to do           */
} twoxyp;

typedef struct Action {
  double     d[6];                  /* Array of doubles describing actions   */
  int        act;                   /* Which action to take?                 */
  int        rcmb;                  /* Recombine spetra after action?        */
  int        nxyp;                  /* Number of pairs of x-y positions      */
  int        nordact;               /* Number of actions made in one go      */
  int        val;                   /* Switch to turn off writing of action  */
  int        i[4];                  /* Array of integers describing actions  */
  twoxyp     *xyp;                  /* Array of pairs of x-y positions       */
} action;

typedef struct CPlot {
  double   swl;                   /* Starting wavelength of large plot       */
  double   ewl;                   /* Ending wavelength of large plot         */
  double   sfac;                  /* Degradation factor for low-res spectrum */
  double   scalclip;              /* Sig-clip level for auto-rescaling orders*/
  double   scalerr;               /* Rel. error threshold for auto-rescaling */
  float    ymx;                   /* Maximum on y-scale                      */
  float    ymn;                   /* Minimum on y-scale                      */
  float    vpd0;                  /* Lower view-port limit for panel 0       */
  float    vpu0;                  /* Upper view-port limit for panel 0       */
  float    vpd1;                  /* Lower view-port limit for panel 1       */
  float    vpu1;                  /* Upper view-port limit for panel 1       */
  float    vpd2;                  /* Lower view-port limit for panel 2       */
  float    vpu2;                  /* Upper view-port limit for panel 2       */
  float    vpd3;                  /* Lower view-port limit for panel 3       */
  float    vpu3;                  /* Upper view-port limit for panel 3       */
  float    vpd4;                  /* Lower view-port limit for panel 4       */
  float    vpu4;                  /* Upper view-port limit for panel 4       */
  int      pgid;                  /* PGPLOT window ID number                 */
  int      sp;                    /* Starting pixel                          */
  int      ep;                    /* Ending pixel                            */
  int      np;                    /* Number of pixels between start and end  */
  int      irank;                 /* Toggle for ranking orders by SNR etc.   */
  int      plval;                 /* Switch for plotting source spectra      */
  int      exit;                  /* Exit code to quit plotting routine      */
  int      refre;                 /* Exit code for refreshing comb plot      */
  int      rccsp;                 /* Exit code for recombination of spectra  */
  int      rscsp;                 /* Exit code for auto-rescaling of spectra */
  char     ptitle[LNGSTRLEN];     /* Plot title                              */
} cplot;

typedef struct RPlot {
  double   zxmin;                 /* Minimum move/zoom factor for x          */
  double   zx;                    /* Move/zoom factor for x                  */
  double   zzx;                   /* Double move/zoom factor for x           */
  double   zymin;                 /* Minimum move/zoom factor for y          */
  double   zy;                    /* Move/zoom factor for y                  */
  double   zzy;                   /* Double move/zoom factor for y           */
  float    ymn;                   /* Min. y for Spec. Nav. panel on this act.*/
  float    ymx;                   /* Max. y for Spec. Nav. panel on this act.*/
  int      sp;                    /* Start pix. for Spec. Nav. on this act.  */
  int      ep;                    /* End pix. for Spec. Nav. on this act.    */
  int      nav;                   /* Does this act. require navig. buttons?  */
  int      actinit;               /* Plot to be initialized for this action? */
  int      cact;                  /* Current action we are working on        */
  int      decx;                  /* Decelerate x flag                       */
  int      decy;                  /* Decelerate y flag                       */
  int      exit;                  /* Exit flag                               */
  int      exun;                  /* Exit/undo toggle                        */
  int      init;                  /* Has the replay plot been initialized?   */
  int      pgid;                  /* PGPLOT window ID of replay plot         */
  int      prre;                  /* Previous/recombine toggle               */
  int      rcmb;                  /* Is a recomb. required after current act?*/
  int      skne;                  /* Skip/next toggle                        */
  int      zomo;                  /* Zoom/move toggle                        */
  char     ptitle[LNGSTRLEN];     /* Plot title                              */
} rplot;

typedef struct Macmap {
  int      nmap;                  /* Number of rows (maps) in Macmap         */
  int      *idx_ori;              /* Array of original indices               */
  int      *idx_new;              /* Array of new indices                    */
  char     mmapfile[NAMELEN];     /* Name of Macmap file                     */
  char     **obj_ori;             /* Array of original names of objects      */
  char     **obj_new;             /* Array of new names of objects           */
  char     **cwl_ori;             /* Array of original central wavelengths   */
  char     **cwl_new;             /* Array of new central wavelengths        */
} macmap;

typedef struct AtMask {
  double   *swl;                  /* Starting wavelengths for mask regions   */
  double   *ewl;                  /* Ending wavelengths for mask regions     */
  int      nmask;                 /* Number of rows in atmosphereic mask file*/
  char     atmaskfile[NAMELEN];   /* Name of mask file                       */
} atmask;

/* LOOKUP TABLES */
/* UVES nominal slit-width--resolving-power product polynomial
   coefficients for the blue and red arms for unbinned and 2x-binned
   settings. These numbers derived from the quality control history
   database at
   http://archive.eso.org/bin/qc1_cgi?action=qc1_browse_table&table=uves_wave
   over from 2010-2016 using slit widths 0.4-1.2" for unbinned
   exposures and 0.8-1.4" for 2x2 binned exposures of the 390 and
   580-nm settings. A simple polynomial fit was performed to derive
   the 2nd order polynomial coefficients. The coefficients produce the
   resolving power, R, at a given slit width (sw) as follows:
   R=a[0]/sw+a[1]+a[2]*sw. There was no discernable difference
   between the redl and redu chips, so the red arm is treated as a
   single entity (values for the redl were used). */
static double UVES_NOM_RES[2][2][3]=
  {{{10033.0,63237.0,-24986.0},  // No binning, blue arm: UVES_NOM_RES[0][0][]
    {8533.3,52709.0,-16005.0}}, // No binning, red arm: UVES_NOM_RES[0][1][]
   {{22011.0,50563.0,-22803.0},  // 2x-binned, blue arm: UVES_NOM_RES[1][0][]
    {28846.0,28505.0,-9533.3}}}; // 2x-binned, red arm: UVES_NOM_RES[1][1][]

/* FUNCTION PROTOTYPES */
double edlen_a2v(double airwl);
double edlen_v2a(double vacwl);
int EW(double *wl, double *fl, double *er, double *co, int np, double *cwl,
       double win, double disp, double res, double *a, double **lim, 
       int *ia, double *ew, double *eew, double *fit, double *chisq,
       int *st, int *niter, int nbtol, int opt, int verbose);
double djmax(double array[], int n);
double djmin(double array[], int n);
float fjmax(float array[], int n);
float fjmin(float array[], int n);
int idxdmax(double *array, int n);
int idxdmin(double *array, int n);
int idxdval(double *array, int n, double val);
int idxfval(float *array, int n, float val);
int UVES_atmask(spectrum *spec, cspectrum *cpsec, atmask *amsk, params *par);
int UVES_boxcar(double *wl, double *fl, double *er, int ndat, int nwid,
		double vwid, double *mean, double *rms, int opt1, int opt2);
int UVES_chunk_cont(double *wl, double *fl, double *er, double *co, int *st,
		    int n, double vc, double rsigl, double rsigu, double pctl,
		    int typ, int ord);
int UVES_combine_cont(spectrum *spec, int nspec, cspectrum *cspec, int opt,
		      params *par);
int UVES_combine_nocont(spectrum *spec, int nspec, cspectrum *cspec, int **comb,
			int ncomb, int csidx, int ceidx, params *par);
int UVES_combine_region(spectrum *spec, int nspec, cspectrum *cspec, params *par);
int UVES_combine_spec(spectrum *spec, int nspec, cspectrum *cspec, params *par);
int UVES_combsynthThAr(cspectrum *cspec, params *par);
int UVES_confit(double *dat, double *err, int *sts, int ndat, int fit_typ,
		int fit_ord, double lrejsig, double urejsig, double pctl,
		int verb, double *fit);
int UVES_cspec_cont(spectrum *spec, int nspec, cspectrum *cspec, params *par);
int UVES_cspec_stats(spectrum *spec, int nspec, cspectrum *cspec, params *par);
int UVES_join_orders(spectrum *spec);
int UVES_hex2dec(char *hex, double *hrs);
int UVES_hirx_blzfit(echorder *ord, double *blz, double initrej, double rejsig,
		     int order, int niter, int opt1, int opt2, double *fit);
int UVES_init_cspec(cspectrum *cspec, params *par, int opt);
int UVES_memspec(spectrum *spec, params *par, int ord, int opt);
int UVES_merge_thar(spectrum *spec, cspectrum *cspec, params *par);
int UVES_model_resol(spectrum *spec);
int UVES_order_cont(spectrum *spec, int nspec, params *par);
int UVES_order_rejsigedge(echorder *ord, params *par);
int UVES_order_sigclip(echorder *ord, params *par);
int UVES_order_stats(echorder *ord, params *par);
int UVES_params_init(params *par);
int UVES_params_set(params *par);
int UVES_past_actions(spectrum *spec, int nspec, cspectrum *cspec, action *act,
		      int nact, int opt, params *par);
double UVES_pixscal(double cwl, double x, int n, int binx);
int UVES_pgenv_init(plotenv *plenv, cplot *cp);
int UVES_plot_cspec(spectrum *spec, int nspec, cspectrum *cspec, cplot *cp,
		    rplot *rp, action **act, int *nact, int *nact_save,
		    params *par);
int UVES_plot_replay(rplot *rp, action *act, int nact, plotbut *but, int nbut,
		     params *par);
int UVES_r1Dspec(cspectrum *cspec, params *par);
int UVES_r2Dspec(spectrum *spec, params *par);
int UVES_r2Dspec_ESOmer(spectrum *spec, params *par);
int UVES_r2Dspec_espresso(spectrum *spec, params *par);
int UVES_r2Dspec_harps(spectrum *spec, params *par);
int UVES_r2Dspec_hirx(spectrum *spec, params *par);
int UVES_r2Dspec_iraf(spectrum *spec, params *par);
int UVES_r2Dspec_iresi(spectrum *spec, params *par);
int UVES_r2Dspec_irls(spectrum *spec, params *par);
int UVES_r2Dspec_KODIAQ(spectrum *spec, params *par);
int UVES_r2Dspec_mage(spectrum *spec, params *par);
int UVES_r2Dspec_makee(spectrum *spec, params *par);
int UVES_ratmask(atmask *amsk, params *par);
int UVES_redispers(spectrum *spec, cspectrum *cspec, long ranseed, params *par);
char *UVES_replace_envinstr(char *str);
int UVES_replay_control(spectrum *spec, int nspec, cspectrum *cspec, cplot *cp,
			rplot *rp, action **act, int *nact, int *nact_save,
			params *par);
int UVES_rescale_region(spectrum *spec, int nspec, cspectrum *cspec, action *act,
			int nact, double scalclip, double scalerr, params *par);
double UVES_revwpol(spectrum *spec, int ord, int idx, double wl, long ranseed,
		    int opt, params *par);
int UVES_rinputfile(char *infile, spectrum **spec, int *nspec, action **act,
		    int *nact, cspectrum *cspec, params *par);
int UVES_rFITSlist(char *infile, spectrum **spec, int *nspec, params *par);
int UVES_rMacmap(macmap *mmap);
int UVES_rscale(spectrum **spec, int nspec, cspectrum *cspec, params *par);
int UVES_rUPLfile(char *infile, spectrum **spec, int *nspec, action **act,
		  int *nact, cspectrum *cspec, params *par);
int UVES_scale(spectrum *spec, int opt);
int UVES_select_subspec(spectrum *spec, int nspec, float wwidth, float wasp,
			params *par);
int UVES_set_wavelen_scale(spectrum *spec, int nspec, cspectrum *cspec,
			   params *par);
int UVES_skysub(spectrum *spec, cspectrum *cspec, params *par);
int UVES_synthThAr(spectrum *spec, params *par);
int UVES_thar_sigarray(spectrum *spec, params *par);
int UVES_undo_lastact(spectrum *spec, int nspec, cspectrum *cspec, action *act,
		      int nact, params *par);
int UVES_vhelio(spectrum *spec);
int UVES_rvshift(spectrum **spec, int nspec, params *par);
int UVES_wDATfile(char *filename, cspectrum *cspec);
int UVES_wFITSfile(spectrum *spec, int nspec, cspectrum *cspec, params *par);
int UVES_wrawFITS(spectrum *spec, params *par);
int UVES_wUPLfile(char *filename, spectrum *spec, int nspec, action *act,
		  int nact, cspectrum *cspec, params *par);
double UVES_wpol(spectrum *spec, int ord, double idx, long ranseed, int opt,
		 params *par);
int UVES_tmpfit(spectrum *spec, params *par);
