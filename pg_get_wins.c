/****************************************************************************
* Subroutine to find ID numbers of all PGPLOT windows
****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "pg_plot.h"
#include "error.h"

int pg_get_wins(Window top, Window pgwins[], int pgwinids[], int *npgwins) {

  int             i;
  unsigned int    nchildren;
  char            *window_name;
  Window          *children, dummy;
  extern Display  *dpy;


  /* Get the name of the top window */
  if (XFetchName(dpy, top, &window_name)) {
    if (strstr(window_name, "PGPLOT Window") != NULL) {
      pgwins[*npgwins] = top;
      if (sscanf(window_name, "%*s %*s %d", &pgwinids[(*npgwins)++]) != 1) {
	nferrormsg("get_pgwins(): Having trouble reading window id from %s!",
		   window_name);
	XFree(window_name);
	return 0;
      }
    }
    else
      XFree(window_name);
  }

  /* Find the children */
  if (!XQueryTree(dpy, top, &dummy, &dummy, &children, &nchildren)) {
    nferrormsg("get_pgwins(): Couldn't query the tree of window top!");
    return 0;
  }

  /* Recurse down the tree */
  for (i = 0; i < nchildren; i++)
    if (!pg_get_wins(children[i], pgwins, pgwinids, npgwins)) {
      XFree(children);
      return 0;
    }

  /* Clean up */
  if (children)
    XFree(children);

  return 1;
}
