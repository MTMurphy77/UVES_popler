/****************************************************************************
* Convert a hexidecimal string (Shh:mm:ss) into decimal representation
****************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "error.h"

int UVES_hex2dec(char *hex, double *hrs) {

  double   hour=0.0,min=0.0,sec=0.0;
  int      sign=-2,nhex=0;

  /* Check first to see if there is a sign at the start (important if
     the input string is a time or RA */
  sign=(isdigit(hex[0])) ? 0 : 1;

  /* Look for a negative sign */
  if (!strncmp(hex,"-",1)) sign=-1;

  /* Read the hours, minutes and seconds from the hexidecimal input */
  if (sign) nhex=sscanf(&(hex[1]),"%lf:%lf:%lf",&hour,&min,&sec);
  else nhex=sscanf(hex,"%lf:%lf:%lf",&hour,&min,&sec);
  if (nhex!=3) {
    nferrormsg("UVES_hex2hrs(): Incorrect format in input string %s.\n\
\tShould have format Shh:mm:ss.s"); return 0;
  }

  /* Convert to hours */
  *hrs=(sign) ? (sec/3600.0+min/60.0+hour)*(double)sign : sec/3600.0+min/60.0+hour;

  return 1;

}
