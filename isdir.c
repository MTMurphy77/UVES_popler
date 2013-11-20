/****************************************************************************
* Check whether file is a directory (not sure how portable this is...)
****************************************************************************/

#include <sys/stat.h>

int isdir(char *name) {

  struct stat    fstat;

  if (stat(name, &fstat)==-1) return 0;
  return S_ISDIR(fstat.st_mode);
}
