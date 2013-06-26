/****************************************************************************
* Bring up a new PGPLOT window which allows the user to select which
* spectra are to be combined. This allows a sub-spectrum to be formed.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "UVES_popler.h"
#include "memory.h"
#include "error.h"

int UVES_select_subspec(spectrum *spec, int nspec, float wwidth, float wasp,
			params *par) {

  float   lstep=0.0,ll=0.0,lu=0.0,lbox=0.0,ltxt=0.0,pgx=0.0,pgy=0.0;
  int     ipage=0,npage=0;
  int     i=0,j=0,k=0,n=0;
  int     *pagen=NULL;
  int     **page=NULL;
  char    pgch[SHORTLEN]="p";
  char    key[LNGSTRLEN]="\0";
  plotenv plenv;
  extern  char   *progname;

  /* Determine required number of pages and allocate memory */
  npage=nspec/NSPPAGE; if (RMDR(nspec,NSPPAGE)) npage++;
  if ((pagen=iarray(npage))==NULL)
    errormsg("UVES_select_subspec(): Cannot allocate memory for pagen\n\
\tarray of size %d",npage);
  if ((page=imatrix(npage,NSPPAGE))==NULL)
    errormsg("UVES_select_subspec(): Cannot allocate memory for page\n\
\tmatrix of size %d",npage,NSPPAGE);

  /* Fill pages with links to spectra */
  for (i=0; i<npage; i++) {
    k=MIN(((i+1)*NSPPAGE),nspec);
    for (j=i*NSPPAGE,pagen[i]=0; j<k; j++) { page[i][pagen[i]]=j; pagen[i]++; }
  }

  /* Initialize plotting */
  UVES_pgenv_init(&plenv,NULL); plenv.wwidth=wwidth; plenv.wasp=wasp; cpgask(0);
  sprintf(key,"%s: Sub-spectrum Selector",progname); pg_open(&plenv,"/XWIN",key,1);
  cpgpage();

  ipage=0; lstep=1.0/(float)NSPPAGE;
  while (strncmp(pgch,"q",1)) {

    /* Begin plotting buffer */
    cpgbbuf(); cpgeras();

    /* Plot main file display window */
    cpgsvp(0.5*plenv.vpl,plenv.vpr,0.7*plenv.vpd,plenv.vpu); cpgsci(0);
    cpgswin(0.0,1.0,0.0,1.0); cpgsci(15); cpgrect(0.0,1.0,0.0,1.0); cpgsci(5);
    cpgslw(3.0*plenv.lw); cpgbox("BC",0.0,0,"BC",0.0,0); cpgsci(3);
    cpgsch(0.7*plenv.ch); cpgslw(plenv.lw); cpgmtxt("T",0.5,0.025,0.5,"Comb?");
    cpgmtxt("T",0.5,0.525,0.5,"Filename"); cpgsci(5);
    cpgmove(0.05,0.0); cpgdraw(0.05,1.0);
    for (i=0,lu=1.0; i<pagen[ipage]; i++) {
      n=ipage*NSPPAGE+i;
      ll=1.0-lstep*(float)(i+1); lbox=0.5*(ll+lu); ltxt=ll+0.4*(lu-ll);
      cpgsci(5); cpgmove(0.0,ll); cpgdraw(1.0,ll);
      if (spec[n].comb) cpgsci(3);
      else cpgsci(2);
      cpgsch(1.5*plenv.ch); cpgpt1(0.025,lbox,-4); cpgsch(0.45*plenv.ch);
      cpgsci(1); cpgptxt(0.06,ltxt,0.0,0.0,spec[n].file); lu=ll;
    }

    /* Plot navigation buttons and print page number */
    cpgsch(1.3*plenv.ch);
    cpgsci(0); cpgsfs(1); cpgslw(plenv.lw); cpgrect(0.0,0.24,-0.12,-0.04); 
    cpgsci(5); cpgsfs(2); cpgslw(2.0*plenv.lw); cpgrect(0.0,0.24,-0.12,-0.04); 
    cpgslw(2.0*plenv.lw); cpgptxt(0.12,-0.06,0.0,0.5,"PREV");
    cpgsci(0); cpgsfs(1); cpgslw(plenv.lw); cpgrect(0.26,0.50,-0.12,-0.04); 
    cpgsci(5); cpgsfs(2); cpgslw(2.0*plenv.lw); cpgrect(0.26,0.50,-0.12,-0.04); 
    cpgslw(2.0*plenv.lw); cpgptxt(0.38,-0.06,0.0,0.5,"NEXT");
    cpgsci(0); cpgsfs(1); cpgslw(plenv.lw); cpgrect(0.75,1.0,-0.12,-0.04); 
    cpgsci(5); cpgsfs(2); cpgslw(2.0*plenv.lw); cpgrect(0.75,1.0,-0.12,-0.04); 
    cpgslw(2.0*plenv.lw); cpgptxt(0.875,-0.06,0.0,0.5,"QUIT");
    sprintf(key,"Page %d of %d",ipage+1,npage); cpgsci(7);
    cpgslw(1.5*plenv.lw); cpgsch(1.0*plenv.ch); cpgptxt(0.625,-0.06,0.0,0.5,key);
    cpgsfs(1);

    /* Flush plotting buffer */
    cpgebuf();

    /* Allow user to make selections */
    cpgsci(7); cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
    if (!strncmp(pgch,"?",1)) {
      fprintf(stderr,"\
Options for sub-spectra selection:\n\
 Left mouse on spectrum    : Toggles selection of that spectrum for combination\n\
 Left mouse on \"Comb?\"   : Toggles selection of all spectra on this page\n\
 a on \"Comb?\"            : Toggles selection of all spectra (on all pages)\n\
 Left mouse on \"PREV\"      : Moves to previous page of spectra\n\
 Left mouse on \"NEXT\"      : Moves to next page of spectra\n\
 q or left mouse on \"QUIT\" : Quits selection and combines selected spectra\n");
    } else if (!strncmp(pgch,"a",1) &&
	       pgx>=0.0 && pgx<=0.05 && pgy>=1.0 && pgy<=1.05) {
      /* Toggle all the combination flags (on all pages) */
      for (n=0; n<nspec; n++) spec[n].comb=(spec[n].comb) ? 0 : 1;
    } else if (!strncmp(pgch,"A",1)) {
      if (pgx>=0.0 && pgx<=1.0 && pgy>=1.0-lstep*(float)pagen[ipage] && pgy<=1.0) {
	/* Determine which spectrum has been clicked on */
	n=ipage*NSPPAGE+(int)((1.0-pgy)/lstep);
	/* Toggle the combination flag for that spectrum */
	spec[n].comb=(spec[n].comb) ? 0 : 1;
      }
      else if (pgx>=0.0 && pgx<=0.05 && pgy>=1.0 && pgy<=1.05) {
	/* Toggle all the combination flags on this page */
	for (i=0; i<pagen[ipage]; i++) {
	  n=ipage*NSPPAGE+i; spec[n].comb=(spec[n].comb) ? 0 : 1;
	}
      }
      else if (pgx>=0.0 && pgx<=0.24 && pgy>=-0.15 && pgy<=-0.01)
	/* Go to previous page */
	ipage=(ipage) ? ipage-1 : npage-1;
      else if (pgx>=0.26 && pgx<=0.50 && pgy>=-0.15 && pgy<=-0.01)
	/* Go to next page */
	ipage=(PMOD(ipage,1,(npage-1),0));
      else if (pgx>=0.75 && pgx<=1.0 && pgy>=-0.15 && pgy<=-0.01)
	/* Quit */
	sprintf(pgch,"%s","q");
    }
  }

  cpgclos();

  /* Clean up */
  free(pagen); free(*page); free(page);

  return 1;

}
