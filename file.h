/***************************************************************************

FILE.H: Include file containing the function prototypes for reading and
writing files.

***************************************************************************/

/* INCLUDE FILES */
#include <stdio.h>

/* DEFINITIONS */

/* STRUCTURES */

/* Prototypes */
FILE    *faskropen(char *query, char *filename, int opt);
FILE    *faskwopen(char *query, char *filename, int opt);
int     isdir(char *);
