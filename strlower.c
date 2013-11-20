/****************************************************************************
Function to convert the entire input string to all lower case letters.
****************************************************************************/

#include <stdio.h>
#include <ctype.h>

int strlower(char *str) {

  char *cptr=NULL;

  for (cptr=str; *cptr; cptr++) *cptr=tolower(*cptr);

  return 1;

}
