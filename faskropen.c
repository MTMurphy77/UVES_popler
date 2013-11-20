/****************************************************************************
* Open a file for reading with optional file name query, checking whether
* the file exists and system escape.
* Options:
* 1 = don't have a default filename
* 2 = have a default, but query the user
* 3 = have a default, query only if the file doesn't exist
* 4 = have a default, return NULL if file doesn't exist
* 5 = have a default, exit program with error message if file doesn't exist
****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "file.h"
#include "input.h"
#include "error.h"

FILE *faskropen(char *query, char *filename, int opt) {

  int status=0;

  if (opt == 1)
    get_input(query, "!%s", filename);
  else if (opt == 2)
    get_input(query, "%s", filename);
  else if (opt == 3 || opt == 4 || opt == 5)
    ;
  else
    errormsg("faskropen(): Unknown option: %d", opt);

  while ((access(filename, R_OK) == -1 || filename[0] == '!') && opt!=4) {
    if (opt == 5)
      errormsg("faskropen(): Cannot open file %s for reading",filename);
    if (filename[0] == '!')
      status=system(filename+1);
    else
      nferrormsg("faskropen(): Cannot open file file %s for reading", filename);
    get_input(query, "!%s", filename);
  }
  return fopen(filename, "r");
}
