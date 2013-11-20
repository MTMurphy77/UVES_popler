/******************************************************************************
* Rename the PGPLOT window with the highest ID
******************************************************************************/

#include <stdio.h>
#include "pg_plot.h"
#include "error.h"

#define N_PGWINS 20 /* Maximum number of PGPLOT windows that can be handled */

Display   *dpy;


int pg_win_rename(char *name) {

  int            screen, pgwinids[N_PGWINS], npgwins;
  char           *display_name=NULL;
  Window         pgwins[N_PGWINS];


  /* Connect to X server */
  if ((dpy = XOpenDisplay(display_name)) == NULL ) {
    nferrormsg("cpgwinrename(): Couldn't connect to X server %s!", 
	       XDisplayName(display_name));
    return 0;
  }
  screen = DefaultScreen(dpy);

  npgwins = 0;
  if (!pg_get_wins(RootWindow(dpy, screen), pgwins, pgwinids, &npgwins)) {
    nferrormsg("cpgwinrename(): Couldn't get PGPLOT window IDs!");
    return 0;
  } else if (npgwins > N_PGWINS) {
    nferrormsg("cpgwinrename(): Too many PGPLOT windows!");
    return 0;
  }

  /* Rename */
  if (npgwins)
    XStoreName(dpy, pgwins[idxival(pgwinids, npgwins,
				   ijmax(pgwinids, npgwins))], name);

  /* Disconnect from X server */
  XCloseDisplay(dpy);

  return 1;
}
