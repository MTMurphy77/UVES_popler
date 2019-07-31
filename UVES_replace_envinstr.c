/****************************************************************************
Function to replace any leading environment variable in a string with
the value for that environment variable
****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

char *UVES_replace_envinstr(char *str) {

  int  l=0,n=0,m=0;
  char *newstr=NULL,*cptr=NULL,*envstr=NULL,*envptr=NULL;

  /* If the input string starts with an environment variable, then try
     to read that variable and replace it with the value. Environment
     variables can be either "${VARIABLE_NAME}" or simply
     "$VARIABLE_NAME" */
  if (!strncmp(str,"$",1) && strlen(str)>=2) {
    /* Test whether variable is enclosed in braces */
    if (!strncmp(&str[1],"{",1)) {
      if (strlen(str)==2) {
	nferrormsg("UVES_replace_envinstr(): Input string appears to start\n\
\twith an environment variable with an opening brace ('{') but is only\n\
\t2 characters long so cannot contain a closing brace ('}'):\n\t%s",str);
	return NULL;
      }
      else {
	cptr=&str[2]; n=strchr(cptr,'}')-cptr; m=1;
      }
    } else {
      cptr=&str[1]; n=strchr(cptr,'/')-cptr; m=0;
    }
    envstr=strndup(cptr,n);
    cptr+=n+m;
    /* Attempt to read the environment variable */
    if ((envptr=getenv(envstr))==NULL) {
      nferrormsg("UVES_replace_envinstr(): Cannot get value of\n\
\tenvironment variable '%s' referred to in input string\n\t'%s'",envstr,str);
      return NULL;
    }
    /* Allocate memory for new string */
    l=strlen(envptr)+strlen(cptr);
    if ((newstr=carray(l))==NULL) {
      nferrormsg("UVES_replace_envinstr(): Cannot allocate memory to string\n\
\tof length %d",l);
      return NULL;
    }
    /* Construct new string */
    strcpy(newstr,envptr); strcat(newstr,cptr);
  } else {
    /* If the input string doesn't start with an environment variable,
       just set the returned pointer to the start of that string */
    newstr=&str[0];
  }

  return newstr;

}
