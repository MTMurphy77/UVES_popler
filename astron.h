/***************************************************************************
ASTRON.H: Include file for operations concerning astronomy
***************************************************************************/

/* INCLUDE FILES */
#include "const.h"

/* DEFINITIONS */
#define AST_MOD(a,b) (b==0.0) ? b : a-((int)(a/b))*b
#define AST_NINT(a) (a-(int)(a)<0.5) ? (int)(a) : (int)(a+1.0)

/* PARAMETERS */
#define AST_RAVSUN    18.0  /* RA of Sun's velocity */
#define AST_DECVSUN   30.0  /* DEC of Sun's velocity */
#define AST_VSUN      20.0  /* Velocity of Sun with respect to LSR (km/s) */
#define AST_EPVSUN  1900.0  /* Epoch of Sun's velocity */

/* STRUCTURES */
typedef struct AstCoord {
  double   epoch;          /* Epoch of coordinates             */
  int      rah;            /* RA hour                          */
  int      ram;            /* RA minute                        */
  double   ras;            /* RA seconds                       */
  double   radec;          /* RA in hours in decimal notation  */
  int      decsign;        /* DEC sign                         */
  int      decd;           /* DEC degree                       */
  int      decm;           /* DEC minute                       */
  double   decs;           /* DEC seconds                      */
  double   decdec;         /* DEC in degrees decimal notation  */
} astcoord;

typedef struct AstObj {
  char     *name;          /* Object name                      */
  astcoord eqbco;          /* Object equiorial B1950 coords    */
  astcoord eqjco;          /* Object equiorial J2000 coords    */
  astcoord galco;          /* Object galactic coords (1950)    */
} astobj;

/* PROTOTYPES */
int ast_coord(double ao, double bo, double ap, double bp, double a1, double b1,
	      double *a2, double *b2);
double ast_date2epoch(int year, int month, int day, double ut);
double ast_date2jd(int year, int month, int day, double t);
int ast_dec2hex(astcoord *coord);
double ast_epoch2jd(double epoch);
int ast_eq2gal(double ra, double dec, double epoch, double *l, double *b);
int ast_hex2dec(astcoord *coord);
int ast_jd2date(double j, int *year, int *month, int *day, double *t);
double ast_jd2epoch(double jd);
double ast_jd2hjd(double ra, double dec, double jd);
double ast_mst(double epoch, double lon);
int ast_precess(double ra1, double dec1, double epoch1, double epoch2,
		double *ra2, double *dec2);
int ast_rotmatrix(double epoch, double **p);
double ast_vbary(double ra, double dec, double epoch);
double ast_vorbit(double ra, double dec, double epoch);
double ast_vr(double ra1, double dec1, double v1, double ra2, double dec2);
double ast_vrotate(double ra, double dec, double epoch, double lat, double lon,
		   double alt);
double ast_vsun(double ra, double dec, double epoch);
