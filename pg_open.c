/****************************************************************************
* Open a graphics device and set plotting environment
* If opt=0 then the default PGPLOT window width and aspect ratio is
* used and recorded in the plotenv structure, if opt=1 then the values
* already within the plotenv structure are used
****************************************************************************/

#include <string.h>
#include "pg_plot.h"
#include "error.h"

int pg_open(plotenv *plenv, char *device, char *progname, int opt) {

  int             i;
  float           cr, cg, cb;
  float           fdum=0.0;
  char            buffer[LNGSTRLEN];

  /* Open the device */
  if (cpgopen(device)<=0)
    errormsg("pg_open(): Couldn't open the device for PGPLOT!"," ");
  cpgask(0);

  /* Enquire colour representation */
  cpgqcr(0,&cr,&cg,&cb);

  /* Is the device interactive? */
  i=LNGSTRLEN; cpgqinf("cursor",buffer,&i);
  if (!strncmp(buffer,"YES",3)) {

    /* Rename the window */
    pg_win_rename(progname);

    /* Use default window size if opt=0 and record that size in
       plotenv structure */
    if (!opt) {
      cpgqvsz(1,&fdum,&(plenv->wwidth),&fdum,&(plenv->wasp));
      plenv->wasp/=plenv->wwidth;
    } else {
      /* Otherwise resize the window */
      cpgpap(plenv->wwidth,plenv->wasp);
    }

    /* Set colours */
    cpgscr(15,cr,cg,cb);
    cpgscr(0,0.184,0.31,0.31);
    plenv->bla=1; plenv->red=2; plenv->gre=3; plenv->blu=4; plenv->cya=5;
    plenv->mag=6; plenv->yel=7;
  }
  else {
    /* Set colours */
    cpgscr(15,cr,cg,cb);
    plenv->bla=1; plenv->red=2; plenv->gre=1; plenv->blu=4; plenv->cya=1;
    plenv->mag=1; plenv->yel=1;
  }

  /* Set other attributes */
  cpgsch(plenv->ch);
  cpgslw(plenv->lw);

  return 1;
}
