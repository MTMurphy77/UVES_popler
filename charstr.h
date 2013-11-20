/***************************************************************************
CHARSTR.H: Include file containing the lengths of characters strings of
           different "types".
***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */
#define TINYLEN          2
#define SHORTLEN         4
#define QUERYLEN         8
#define NAMELEN         64
#define LNGSTRLEN      128
#define VLNGSTRLEN     256
#define VVLNGSTRLEN    512
#define VVVLNGSTRLEN  1024
#define HUGESTRLEN    2048
#define VHUGESTRLEN   4096
#define VVHUGESTRLEN  8192
#define VVVHUGESTRLEN 16384
#define UBERSTRLEN    32768

/* STRUCTURES */

/* FUNCTION PROTOTYPES */
int strisnum(char *string, int length, int opt);
int strlower(char *str);
int strupper(char *str);
char *strtolower(char *string, int length);
char *strtoupper(char *string, int length);
