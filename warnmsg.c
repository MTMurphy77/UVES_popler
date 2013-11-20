/***************************************************************************
* Send a warning message to the terminal
***************************************************************************/

#include <stdio.h>
#include <stdarg.h>

void warnmsg(char *fmt, ...) {

  va_list       args;
  extern char   *progname;

  va_start(args, fmt);
  fprintf(stderr, "%s: WARNING: ", progname);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");
}
