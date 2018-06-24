/****************************************************************************
* Get a string from stdin, character by character, with file name completion
* and history; a maximum of size characters are read, no trailing newline.
****************************************************************************/

#include <termios.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "charstr.h"
#include "input.h"
#include "error.h"

/****************************************************************************
* A function to set some terminal properties;
* Hacked from http://www.ohse.de/uwe/articles/getch.html
****************************************************************************/

static int io_tio_set_flag_c_lflag(int fd, int val, int on, int *old) {
  
  struct termios tio;

  if (tcgetattr(fd, &tio))
    return 1;
  if (old) 
    *old = (tio.c_lflag & (val));
  if (on)
    tio.c_lflag |= (val);
  else
    tio.c_lflag &= ~(val);
  if (tcsetattr(fd, TCSADRAIN, &tio))
    return 1;
  return 0;
}

/*
static int io_tio_get_flag_##which(int fd, int bit, int *value) \
{ struct termios tio; \
    if (tcgetattr(fd,&tio)) return -1; \
        *value=(tio.which & (bit)); \
    return 0; \
}
*/


/***************************************************************************/


int getscbc(char *s, int size) {

  int          i, old_icanon, old_echo;
  static int   ihist=0, nhist=0;
  char         *buf;
  static char  hist[N_HIST][VLNGSTRLEN];

  /* Turn off line by line mode */
  if (io_tio_set_flag_c_lflag(0, ICANON, 0, &old_icanon)) {
    nferrormsg("getscbc(): Couldn't turn off ICANON!");
    return 1;
  }

  /* Turn off echoing */
  if (io_tio_set_flag_c_lflag(0, ECHO, 0, &old_echo)) {
    nferrormsg("getscbc(): Couldn't turn off ECHO!");
    return 1;
  }

  i = 0;
  while (i < size && read(0, &s[i], 1) > 0) {
    if (s[i] == '\n' || s[i] == '\r' || s[i]+64 == 'D') {
      /* End of input */
      break;

    } else if (s[i] == 127) {
      /* Backspace */
      if (i > 0) {
	fprintf(stderr, "\b \b");
	i--;
      }
      continue;

    } else if (s[i] == '\t') {
      /* File name completion */
      s[i] = '\0';
      buf = ((buf = strrchr(s, ' ')) != NULL) ? buf+1 : s;
      fcompl(buf);
      i = strlen(s);
      continue;

    } else if (s[i] == '\033') {
      /* Escape sequence */
      read(0, &s[i], 1);
      if (s[i] == '[') {
	read(0, &s[i], 1);

	if (s[i] == 'A') {
	  /* History */
	  if (ihist == 0)
	    fprintf(stderr, "\a");
	  else {
	    s[i] = '\0';
	    if (ihist == nhist) 
	      strcpy(hist[ihist], s);
	    if ((i = strlen(s)) != 0) {
	      sprintf(s, "\033[%dD\033[%dP", i, i);
	      fprintf(stderr, "%s", s);
	    }
	    strcpy(s, hist[--ihist]);
	    fprintf(stderr, "%s", s);
	    i = strlen(s);
	  }
	} else if (s[i] == 'B') {
	  /* History */
	  if (ihist == nhist)
	    fprintf(stderr, "\a");
	  else {
	    s[i] = '\0';
	    if ((i = strlen(s)) != 0) {
	      sprintf(s, "\033[%dD\033[%dP", i, i);
	      fprintf(stderr, "%s", s);
	    }
	    strcpy(s, hist[++ihist]);
	    fprintf(stderr, "%s", s);
	    i = strlen(s);
	  }
	}
      }
      continue;

    } else if (s[i] < 32) {
      /* Other control characters */
      continue;

    } else {
      fprintf(stderr, "%c", s[i]);
    }
    i++;
  }
  fprintf(stderr, "\n");
  s[i] = '\0';

  /* Push on history stack */
  strcpy(hist[nhist], s);
  if (nhist != N_HIST-1)
    nhist++;  
  else {
    for (i = 1; i < N_HIST; i++)
      strcpy(hist[i-1], hist[i]);
  }
  ihist = nhist;

  /* Restore input and echo modes */
  io_tio_set_flag_c_lflag(0, ICANON, old_icanon, NULL);
  io_tio_set_flag_c_lflag(0, ECHO, old_echo, NULL);

  return 0;
}
