/***************************************************************************
* Send a fatal error message to the terminal and exit
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void errormsg(char *fmt, ...) {

  va_list       args;
  extern char   *progname;

  va_start(args, fmt);
  fprintf(stderr, "\a%s: FATAL ERROR: ", progname);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");
  exit(1);
}
