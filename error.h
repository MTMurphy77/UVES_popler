/***************************************************************************

ERROR.H: Include file containing the function prototypes handling errors

***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */
#define MAXNERR 3  /* Maximum number of errors allowed */

/* STRUCTURES */

/* Prototypes */
void    errormsg(char *fmt, ...);
void    nferrormsg(char *fmt, ...);
void    warnmsg(char *fmt, ...);
