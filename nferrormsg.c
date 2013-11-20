/***************************************************************************
* Send a non-fatal error message to the terminal
***************************************************************************/

#include <stdio.h>
#include <stdarg.h>

void nferrormsg(char *fmt, ...) {

  va_list       args;
  extern char   *progname;

  va_start(args, fmt);
  fprintf(stderr, "\a%s: ERROR: ", progname);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");
}
