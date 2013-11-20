/****************************************************************************
* Complete a file name; returns 0 on successful (partial) completion
****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include "charstr.h"
#include "file.h"
#include "error.h"

void fcompl(char *oname) {

  int            i, n, nmatch;
  char           *buf, name[VLNGSTRLEN], dir[VLNGSTRLEN], match[VLNGSTRLEN];
  DIR            *dirp;
  struct dirent  *dp;

  /* Get the directory and (incomplete) file name */
  strcpy(dir, oname);
  if ((buf = strrchr(dir, '/')) != NULL) {
    strcpy(name, buf+1);
    *buf = '\0';
  } else {
    strcpy(name, dir);
    strcpy(dir, ".");
  }

  /* Open the directory */
  if ((dirp = opendir(dir)) == NULL) {
    nferrormsg("fcompl(): Couldn't open directory %s!", name);
    return;
  }

  /* Find matches */
  n = strlen(name);
  nmatch = 0;
  while ((dp = readdir(dirp)) != NULL)
    if (!strncmp(dp->d_name, name, n)) {
      if (nmatch == 0) {
	/* First match */
	nmatch = 1;
	strcpy(match, dp->d_name);
      } else {
	/* More than one match */
	i = n + 1;
	while (!strncmp(dp->d_name, match, i))
	  i++;
	match[i-1] = '\0';
	nmatch++;
      }
    }
  closedir(dirp);

  if (nmatch == 0) {
    fprintf(stderr, "\a");
    return;
  }

  if (nmatch > 1)
    fprintf(stderr, "\a");
  else {
    /* If there was only one match, is it a directory? */
    strcat(dir, "/");
    if (isdir(strcat(dir, match)))
      strcat(match, "/");
    else
      strcat(match, " ");
  }

  strcat(oname, match+n);
  fprintf(stderr, "%s", match+n);
}
