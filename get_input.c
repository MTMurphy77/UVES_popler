/****************************************************************************
* Get input from STDIN.
* Only integers, floats, doubles, characters and strings are allowed.
* Values stored in the arguments will be offered as defaults unless
* the first character of the format string is not a '%'.
****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "charstr.h"
#include "error.h"
#include "input.h"

void get_input(char *query, char *fmt, ...) {

  int         *ip, i, flag, error=0;
  float       *fp;
  double      *dp;
  char        *p, *cp, deffmt[VLNGSTRLEN], answer[VLNGSTRLEN], *pans;
  va_list     ap;
  
  /* Initialise argument list */
  va_start(ap, fmt);

  /* Print the query */
  p = fmt;
  /* Check for defaults */
  if (*p != '%')
    fprintf(stderr, "%s ", query);
  else {
    fprintf(stderr, "%s [", query);
    
    /* Step through the format string */
    for (; *p; p++) {

      /* Find '%', print everything else */
      if (*p != '%') {
	fputc(*p, stderr);
	continue;
      }

      /* Get extra formatting information if present */
      p++;
      i = 0;
      deffmt[i++] = '%';
      while (*p == '0' || *p == '1' || *p == '2' || *p == '3' || *p == '4' ||
	     *p == '5' || *p == '6' || *p == '7' || *p == '8' || *p == '9' || 
	     *p == '.')
	deffmt[i++] = *p++;

      /* Print defaults */
      switch (*p) {
      case 'd':
	deffmt[i++] = 'd';
	deffmt[i] = '\0';
	ip = va_arg(ap, int *);
	fprintf(stderr, deffmt, *ip);
	break;
      case 'f':
	deffmt[i++] = 'f';
	deffmt[i] = '\0';
	fp = va_arg(ap, float *);
	fprintf(stderr, deffmt, *fp);
	break;
      case 'e':
	deffmt[i++] = 'e';
	deffmt[i] = '\0';
	fp = va_arg(ap, float *);
	fprintf(stderr, deffmt, *fp);
	break;
      case 'g':
	deffmt[i++] = 'g';
	deffmt[i] = '\0';
	fp = va_arg(ap, float *);
	fprintf(stderr, deffmt, *fp);
	break;
      case 'l':
	deffmt[i++] = *++p;
	deffmt[i] = '\0';
	/*
	p++;
	*/
	dp = va_arg(ap, double *);
	fprintf(stderr, deffmt, *dp);
	break;
      case 'c':
	deffmt[i++] = 'c';
	deffmt[i] = '\0';
	cp = va_arg(ap, char *);
	fprintf(stderr, deffmt, *cp);
	break;
      case 's':
	deffmt[i++] = 's';
	deffmt[i] = '\0';
	cp = va_arg(ap, char *);
	fprintf(stderr, deffmt, cp);
	break;
      default:
	errormsg("get_input(): Don't know format type!");
	break;
      }
    }
    fprintf(stderr, "]: ");
  }
  
  flag = 1;
  while (flag) {
    
    /* Scan the data */
    getscbc(answer, VLNGSTRLEN);

    if (strlen(answer) == 0 && fmt[0] == '%') {
      /* Defaults are present and were selected; clean up and leave */
      va_end(ap);
      return;
    }
    
    /* Initialise again */
    va_end(ap);
    va_start(ap, fmt);
      
    /* Step through format string */
    flag = 0;
    for (p = fmt, pans = answer; *p; p++) {

      /* Find '%', skip everything else */
      if (*p != '%')
	continue;
      
      /* Get extra formatting information if present (only used for
	 characters) */
      p++;
      i = 0;
      deffmt[i++] = '%';
      while (*p == '0' || *p == '1' || *p == '2' || *p == '3' || *p == '4' ||
	     *p == '5' || *p == '6' || *p == '7' || *p == '8' || *p == '9' ||
	     *p == '.')
	deffmt[i++] = *p++;
	
      /* Read in data, flag errors */
      switch (*p) {
      case 'd':
	ip = va_arg(ap, int *);
	if (sscanf(pans, "%d", ip) != 1)
	  flag = 1;
	break;
      case 'f':
	fp = va_arg(ap, float *);
	if (sscanf(pans, "%f", fp) != 1)
	  flag = 1;
	break;
      case 'e':
	fp = va_arg(ap, float *);
	if (sscanf(pans, "%f", fp) != 1)
	  flag = 1;
	break;
      case 'g':
	fp = va_arg(ap, float *);
	if (sscanf(pans, "%f", fp) != 1)
	  flag = 1;
	break;  
      case 'l':
	p++;
	dp = va_arg(ap, double *);
	if (sscanf(pans, "%lf", dp) != 1)
	  flag = 1;
	break;
      case 'c':
	deffmt[i++] = 'c';
	deffmt[i] = '\0';
	cp = va_arg(ap, char *);
	if (sscanf(pans, deffmt, cp) != 1)
	  flag = 1;
	break;
      case 's':
	cp = va_arg(ap, char *);
	if (sscanf(pans, "%s", cp) != 1)
	  flag = 1;
	break;
      default:
	errormsg("get_input(): Don't know format type!");
	break;
      }

      /* Move to next answer */
      pans = strchr(pans, ' ');
      if (pans == 0) {
	if (fmt[0] != '%' && *(p+1) != '\0')
	  flag = 1;
	break;
      } else 
	pans++;
    }

    /* Repeat input if error occured */
    if (flag && error++<MAXNERR) {
      fprintf(stderr, "get_input(): Bad input, try again!\n");
      p=fgets(answer, VLNGSTRLEN, stdin);
    }
    /* Exit with error message if too many errors have occurred */
    else if (error>=MAXNERR)
      errormsg("get_input(): Too many errors!");
  }

  /* Clean up */
  va_end(ap);
}
