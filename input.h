/***************************************************************************

INPUT.H: Include file for using subroutines like get_input() etc.

***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */
#define N_HIST   11   /* Number of lines kept in input history */

/* STRUCTURES */

/* Prototypes */
void get_input(char *query, char *fmt, ...);
int  getscbc(char *, int);
void fcompl(char *);
