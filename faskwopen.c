/****************************************************************************
* Open a file for writing with optional file name query, checking whether
* the file exists and system escape.
* Options:
* 1 = don't have a default
* 2 = have a default but query the user
* 3 = have a default, query only if the file exists
* 4 = have a default, use it even if the file exists (like fopen)
* 5 = have a default, exit program with error message if it doesn't exist
* 6 = as 4, but append to file instead of overwriting file
****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "file.h"
#include "input.h"
#include "error.h"

FILE *faskwopen(char *query, char *filename, int opt) {

  char    ans;
  
  if (opt == 1 || opt == 2 || opt == 3) {

    if (opt == 1) get_input(query, "!%s", filename);
    else if (opt == 2) get_input(query, "%s", filename);

    while (!access(filename,W_OK) || filename[0] == '!') {
      if (filename[0] == '!') system(filename+1);
      else {
	warnmsg("File %s exists!", filename);
	ans = 'y';
	get_input("Overwrite (y/n)?", "%c", &ans);
	if (ans == 'y')
	  break;
      }
      get_input(query, "!%s", filename);
    }
  }
  else if (opt == 4)
    ;
  else if (opt == 5 && access(filename,W_OK))
    errormsg("faskwopen(): Cannot open file %s for writing",filename);
  else if (opt == 6) return fopen(filename, "a");
  else errormsg("faskwopen(): Unknown option: %d", opt);

  return fopen(filename, "w");

}
