/************************************************************************** 
* Draw a button with PGPLOT
**************************************************************************/

#include "pg_plot.h"

int pg_button(plotbut *but) {

  int   oi=0,ofs=0;
  float och=0.0,xch=0.0,ych=0.0;

  /* Get current values of colour index, character height and fill style */
  cpgqci(&oi); cpgqch(&och); cpgqfs(&ofs);

  /* Plot button */
  cpgsci(but->bgc); cpgsfs(1);
  cpgrect((float)but->x1,(float)but->x2,(float)but->y1,(float)but->y2);
  cpgsci(but->boc); cpgsfs(2);
  cpgrect((float)but->x1,(float)but->x2,(float)but->y1,(float)but->y2);
  cpgsci(but->fgc); cpgsch(but->ch); cpgqcs(4,&xch,&ych);
  cpgptxt(0.5*(float)(but->x1+but->x2),0.5*(float)(but->y1+but->y2)-0.35*ych,0.0,
	  0.5,but->lab);

  /* Reset colour index, character height and fill style */
  cpgsci(oi); cpgsch(och); cpgsfs(ofs);

  return 1;

}
