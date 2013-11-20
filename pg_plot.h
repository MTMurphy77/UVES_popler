/***************************************************************************

PG_PLOT.H: Include file containing information relevant to producing plots
in PGPLOT

***************************************************************************/

/* INCLUDE FILES */
#include <cpgplot.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/Xos.h>
#include <X11/extensions/shape.h>
#include "charstr.h"
#include "utils.h"

/* DEFINITIONS */
#define   PGNSUB      12
#define   PGBLA        1
#define   PGRED        2
#define   PGGRE        3
#define   PGBLU        4
#define   PGCYA        5
#define   PGMAG        6
#define   PGYEL        7

#define   POINT        0
#define   HLINE        1
#define   VLINE        2

/* STRUCTURES */
typedef struct PlotEnv {
  float    wwidth;                   /* Window width */
  float    wasp;                     /* Window aspect ratio */
  float    ch;                       /* Character height */
  float    vpu;                      /* View port upper limit */
  float    vpd;                      /* View port lower limit */
  float    vpl;                      /* View port left limit */
  float    vpr;                      /* View port right limit */
  float    xmin[PGNSUB];             /* Minimum x-value for each subpanel */
  float    xmax[PGNSUB];             /* Maximum x-value for each subpanel */
  float    ymin[PGNSUB];             /* Minimum y-value for each subpanel */
  float    ymax[PGNSUB];             /* Maximum y-value for each subpanel */
  int      lw;                       /* Line width */
  int      nxsub;                    /* Number of subpanels in x-direction */
  int      nysub;                    /* Number of subpanels in y-direction */
  int      bla;                      /* Index to be used for BLACK         */
  int      red;                      /* Index to be used for RED           */
  int      gre;                      /* Index to be used for GREEN         */
  int      blu;                      /* Index to be used for BLUE          */
  int      cya;                      /* Index to be used for CYAN          */
  int      mag;                      /* Index to be used for MAGENTA       */
  int      yel;                      /* Index to be used for YELLOW        */
  char     xlab[PGNSUB][LNGSTRLEN];  /* Label for x-axis for each subpanel */
  char     ylab[PGNSUB][LNGSTRLEN];  /* Label for y-axis for each subpanel */
  char     title[PGNSUB][LNGSTRLEN]; /* Title for each subpanel */
} plotenv;

typedef struct PlObj {
  double x;                          /* X position of object               */
  double y;                          /* Y position of object               */
  int    m;                          /* Object "morphology", i.e. hline etc*/
  int    s;                          /* PGPLOT symbol used to represent obj*/
  int    c;                          /* PGPLOT colour used for object      */
} plobj;

typedef struct PlotBut {
  double ch;                         /* Character height                   */
  double x1;                         /* x-location of first corner         */
  double x2;                         /* x-location of second corner        */
  double y1;                         /* y-location of first corner         */
  double y2;                         /* y-location of second corner        */
  int    bgc;                        /* Background colour index            */
  int    boc;                        /* Border colour index                */
  int    fgc;                        /* Foreground colour index            */
  char   lab[LNGSTRLEN];             /* Label on button                    */
} plotbut;

/* PROTOTYPES */
int pg_addobj(double x, double y, int m, int s, int c, plobj **obj, int *nobj);
int pg_button(plotbut *but);
int pg_open(plotenv *plenv, char *device, char *progname, int opt);
int pg_palette(int type, float contrast, float brightness);
int pg_palette_fiddle(int *type, float *contrast, float *brightness);
int pg_win_rename(char *name);
int pg_get_nearobj(float xpos, float ypos, double xnorm, double ynorm,
		   plobj *obj, int nobj, int opt);
int pg_get_wins(Window top, Window pgwins[], int pgwinids[], int *npgwins);
