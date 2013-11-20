/****************************************************************************
Function to decide if a string contains only numeric data, i.e. only
digits, decimal points and "+" or "-" signs. If opt=0, scientific notation
(e.g. 1.e-6 or 1.E+6) is allowed. If opt=1, it is not allowed.
****************************************************************************/

#include <ctype.h>
#include "charstr.h"
#include "error.h"

int strisnum(char *string, int length, int opt) {

  int result=0;
  int i=0;

  for (i=0; i<length; i++) {
    if (string[i]=='\0' || string[i]=='\n') break;
    if (!isdigit(string[i]) && string[i]!='+' && string[i]!='-' &&
	string[i]!='.' && string[i]!=' ') {
      if (opt) return result;
      if (string[i]=='e' || string[i]=='E') {
	if (i+1==length || string[i+1]=='\0' || string[i+1]=='\n') return result;
        if (isdigit(string[i]) && string[i]!='+' && string[i]!='-' &&
	    string[i]!='.' && string[i]!=' ') return result;
      } else return result;
    }
  }
  result=1;

  return result;
    
}
