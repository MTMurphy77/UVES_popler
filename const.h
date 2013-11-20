/****************************************************************************
* Constants
****************************************************************************/

#include <math.h>

/* Numbers */
#define C_E              M_E                    /* e */
#define C_LOG2E          M_LOG2E                /* log_2 e */
#define C_LOG10E         M_LOG10E               /* log_10 e */
#define C_LN2            M_LN2                  /* log_e 2 */
#define C_LN10           M_LN10                 /* log_e 10 */
#define C_PI             M_PI                   /* pi */
#define C_PI_2           M_PI_2                 /* pi/2 */
#define C_PI_4           M_PI_4                 /* pi/4 */
#define C_1_PI           M_1_PI                 /* 1/pi */
#define C_2_PI           M_2_PI                 /* 2/pi */
#define C_2_SQRTPI       M_2_SQRTPI             /* 2/sqrt(pi) */
#define C_SQRT2          M_SQRT2                /* sqrt(2) */
#define C_SQRT3          1.73205080756887729352 /* sqrt(3) */
#define C_SQRT5          2.23606797749978969640 /* sqrt(5) */
#define C_1_SQRT2        M_SQRT1_2              /* 1/sqrt(2) */
#define C_2PI            6.28318530717958647692 /* 2*pi */
#define C_4PI           12.5663706143591729538  /* 4*pi */
#define C_SQRTPI         1.772453850905516      /* Square root of pi */
#define C_RPDEG          0.017453292519943      /* Radians/degree: pi/180 */
#define C_AMPDEG        60.0                    /* Arcminutes per degree */
#define C_ASPDEG      3600.0                    /* Arcseconds per degree */

/* Physics */
#define C_C      299792458.0                    /* Speed of light [m/s] */
#define C_C_K       299792.458                  /* Speed of light [km/s] */
#define C_PC             3.085678e16            /* Parsec [m] */
#define C_KPC            3.085678e19            /* Kiloparsec [m] */
#define C_MPC            3.085678e22            /* Megaparsec [m] */
#define C_FSC            0.0072973527442        /* Fine structure constant */
#define C_RFSC         137.03599958             /* Reciprocal FSC */

/* Astronomy and astronomical time */
#define C_SECMIN        60.0                    /* Seconds in a minute */
#define C_MINHR         60.0                    /* Minutes in an hour */
#define C_SECHR       3600.0                    /* Seconds in an hour */
#define C_SECDAY     86400.0                    /* Seconds in a (Julian) day */
#define C_SECYEAR 31557600.0                    /* Seconds in (Julian) year */
#define C_J2000       2000.0                    /* J2000 epoch */
#define C_JD2000   2451545.0                    /* J2000 Julian Date */
#define C_JYEAR        365.25                   /* Julian year */

/* Astronomical coordinate systems */
#define C_LONGNCP      123.00                   /* Longtitude of NCP (1950) */
#define C_RAGPOLE      192.25                   /* RA Galactic Pole (12:49, 1950)*/
#define C_DECGPOLE      27.4                    /* DEC Galactic Pole (27:24,1950)*/
#define C_GEPOCH      1950.0                    /* Epoch of galactic coord defn */

/* Other */
#define C_FWHMSIG 2.354820045030949327          /* Ratio Gaussian FWHM to sigma */
