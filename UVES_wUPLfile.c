/****************************************************************************
* Write out in a UVES_popler log file to record input file list and any
* historical actions
****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "UVES_popler.h"
#include "file.h"
#include "error.h"

int UVES_wUPLfile(char *filename, spectrum *spec, int nspec, action *act,
		  int nact, cspectrum *cspec, params *par) {

  double    version=0.0;
  int       nval=0;
  int       i=0,j=0;
  char      user[NAMELEN],tstr[LNGSTRLEN];
  char      flxerrstr[NAMELEN]="\0",befaftstr[NAMELEN]="\0";
  FILE      *data_file=NULL;
  struct tm *tptr;
  time_t    lt;

  /* Open specified file */
  if (!access(filename,R_OK))
    warnmsg("Overwriting existing file %s",filename);
  else fprintf(stderr,"Opening UVES_popler log file %s\n",filename);
  if ((data_file=faskwopen("UVES_popler log file name?",filename,4))==NULL)
    errormsg("Can not open file %s",filename);

  /* Decide which version number to write to output UPL file, either
     the version number in the input UPL file or the version number of
     the running version of UVES_popler */
  version=(par->backvers) ? par->version : VERSION;

  /* First line of log file is a timestamp and user info */
  sprintf(user,"%s",getenv("USER")); lt=time(NULL); tptr=localtime(&lt);
  // strftime(tstr,LNGSTRLEN,"%a %b %d %Y %H:%M:%S %Z (UTC%z)",tptr);
  strftime(tstr,LNGSTRLEN,"%Y-%m-%d:%H:%M:%S(UTC%z)",tptr);
  fprintf(data_file,"UVES_popler : Version %.2lf UPL (v%.2lf code) : %s : %s\n",version,VERSION,tstr,user);
  fprintf(data_file,"\n");

  /* Print reduction parameters */
  fprintf(data_file,"General spectrum specifications:\n");
  fprintf(data_file," combmeth = %d\n",par->combmeth);
  fprintf(data_file," disp = %17.12lf\n",par->disp);
  fprintf(data_file," filetype = %d\n",par->filetype);
  if (version>=0.31) {
    fprintf(data_file," helio = %d\n",par->helio);
    fprintf(data_file," vacwl = %d\n",par->vacwl);
  }
  fprintf(data_file," zem = %10.7lf\n",par->zem);
  fprintf(data_file," linear = %d\n",par->linear);
  fprintf(data_file," thar = %d\n",par->thar);
  fprintf(data_file,"Initial sigma-clipping of orders:\n");
  fprintf(data_file," nordclip = %d\n",par->nordclip);
  fprintf(data_file," nordsig = %d\n",par->nordsig);
  fprintf(data_file," ordsig = %9.4lf\n",par->ordsig);
  fprintf(data_file," ordsignbr = %9.4lf\n",par->ordsignbr);
  fprintf(data_file," ordsigzero = %9.4lf\n",par->ordsigzero);
  fprintf(data_file," ordmedfrac = %9.4lf\n",par->ordmedfrac);
  fprintf(data_file," ordmedrej = %9.4lf\n",par->ordmedrej);
  fprintf(data_file,"Spectrum/order combination parameters:\n");
  fprintf(data_file," clipsig = %9.4lf\n",par->clipsig);
  fprintf(data_file," lrgerr = %9.4lf\n",par->lrgerr);
  fprintf(data_file," rankspec = %d\n",par->rankspec);
  fprintf(data_file," scalmeth = %d\n",par->scalmeth);
  fprintf(data_file," scalclip = %9.4lf\n",par->scalclip);
  fprintf(data_file," scalerr = %9.4lf\n",par->scalerr);
  fprintf(data_file," nscalclip = %d\n",par->nscalclip);
  fprintf(data_file,"Order continuum fitting method parameters:\n");
  fprintf(data_file," contftyp = %d\n",par->contftyp);
  fprintf(data_file," contord = %d\n",par->contord);
  fprintf(data_file," contpctl = %9.4lf\n",par->contpctl);
  fprintf(data_file," contsigl = %9.4lf\n",par->contsigl);
  fprintf(data_file," contsigu = %9.4lf\n",par->contsigu);
  fprintf(data_file," contwgt = %d\n",par->contwgt);
  fprintf(data_file,"Combined spectrum continuum fitting method parameters:\n");
  fprintf(data_file," cordlya = %d\n",par->cordlya);
  fprintf(data_file," cordred = %d\n",par->cordred);
  fprintf(data_file," ftyplya = %d\n",par->ftyplya);
  fprintf(data_file," ftypred = %d\n",par->ftypred);
  fprintf(data_file," nocont = %d\n",par->nocont);
  fprintf(data_file," pctllya = %9.4lf\n",par->pctllya);
  fprintf(data_file," pctlred = %9.4lf\n",par->pctlred);
  fprintf(data_file," rsiglyal = %9.4lf\n",par->rsiglyal);
  fprintf(data_file," rsiglyau = %9.4lf\n",par->rsiglyau);
  fprintf(data_file," rsigredl = %9.4lf\n",par->rsigredl);
  fprintf(data_file," rsigredu = %9.4lf\n",par->rsigredu);
  fprintf(data_file," vclya = %9.4lf\n",par->vclya);
  fprintf(data_file," vcred = %9.4lf\n",par->vcred);
  fprintf(data_file," vlya = %9.4lf\n",par->vlya);
  fprintf(data_file,"Data file options:\n");
  fprintf(data_file," dat = %d\n",par->dat);
  if (version>0.52) fprintf(data_file," raw = %d\n",par->raw);
  fprintf(data_file," save = %d\n",par->save);
  if (version>=0.56) fprintf(data_file," distort = %d\n",par->distort);
  fprintf(data_file,"\n");

  /* Print number of files and then list each path and absolute path name */
  fprintf(data_file,"Number_of_spectra = %d\n",nspec);
  for (i=0; i<nspec; i++) {
    if (par->filetype==FTMIX) fprintf(data_file,"%3.3d_FTYP = %d\n",i+1,spec[i].ftype);
    fprintf(data_file,"%3.3d_PATH = %s\n",i+1,spec[i].path);
    fprintf(data_file," %3.3d_FLUX = %s\n",i+1,spec[i].abfile);
    if (par->thar<=1) fprintf(data_file," %3.3d_ERRO = %s\n",i+1,spec[i].aberfile);
    if (par->thar==1) fprintf(data_file," %3.3d_THAR = %s\n",i+1,spec[i].abthfile);
    fprintf(data_file," %3.3d_WPOL = %s\n",i+1,spec[i].abwlfile);
    if (version>=0.73) {
      for (j=0; j<spec[i].nor; j++) {
	if (spec[i].or[j].insclbefaft) {
	  if (spec[i].or[j].insclfe==0) sprintf(flxerrstr,"f");
	  else if (spec[i].or[j].insclfe==1) sprintf(flxerrstr,"e");
	  else sprintf(flxerrstr,"fe");
	  if (spec[i].or[j].insclbefaft==1) sprintf(befaftstr,"b");
	  else sprintf(befaftstr,"a");
	  fprintf(data_file," %3.3d_SCAL_%2.2d = %s %lf %s\n",i+1,j+1,flxerrstr,
		  spec[i].or[j].inscl,befaftstr);
	}
      }
    }
    if (version>=0.67) fprintf(data_file," %3.3d_VSHT = %lf %lf %lf\n",i+1,
			       spec[i].vshift,spec[i].vslope,spec[i].refwav);
    else if (version>=0.50) fprintf(data_file," %3.3d_VSHT = %lf\n",i+1,spec[i].vshift);
    if (version>=0.56) fprintf(data_file," %3.3d_DTRT = %ld\n",i+1,spec[i].distort_seed);
  }
  if (par->filetype==FTCOMB) {
    fprintf(data_file,"%3.3d_PATH = %s\n",i,cspec->path);
    fprintf(data_file," %3.3d_COMB = %s\n",i,cspec->abfile);
  }
  fprintf(data_file,"\n");

  /* Determine number of valid actions */
  for (i=0; i<nact; i++) if (act[i].val) nval++;

  /* Write out historical actions */
  fprintf(data_file,"Number_of_actions = %d\n",nval);
  for (i=0; i<nact; i++) {
    if (act[i].val) {
      fprintf(data_file,"%3.3d_ACTN = %d  %3.3d_RCMB = %d\n",i+1,act[i].act,
	      i+1,act[i].rcmb);
      if (act[i].act<=ICACT)
	fprintf(data_file," %3.3d_CORD = %14.8lf %14.8lf %14.8lf %14.8lf\n",
		i+1,act[i].d[0],act[i].d[1],act[i].d[2],act[i].d[3]);
      switch (act[i].act) {
      case COACT:
	fprintf(data_file," %3.3d_SPEC = %d  %3.3d_ORDR = %d\n",
		i+1,act[i].i[0]+1,i+1,act[i].i[1]+1);
	break;
      case UOACT:
	fprintf(data_file," %3.3d_SPEC = %d  %3.3d_ORDR = %d\n",
		i+1,act[i].i[0]+1,i+1,act[i].i[1]+1);
	break;
      case FOACT:
	fprintf(data_file," %3.3d_SPEC = %d  %3.3d_ORDR = %d  %3.3d_FTYP = %d  \
%3.3d_FORD = %d  %3.3d_REJL = %5.2lf  %3.3d_REJU = %5.2lf  %3.3d_NXYP = %d\n",
		i+1,act[i].i[0]+1,i+1,act[i].i[1]+1,i+1,act[i].i[2],i+1,
		act[i].i[3],i+1,act[i].d[4],i+1,act[i].d[5],i+1,act[i].nxyp);
	for (j=0; j<act[i].nxyp; j++) {
	  fprintf(data_file,
		  "  %3.3d_XY_%3.3d = %d  %14.8lf  %14.8lf  %14.8lf  %14.8lf\n",
		  i+1,j+1,act[i].xyp[j].i,act[i].xyp[j].x1,act[i].xyp[j].x2,
		  act[i].xyp[j].y1,act[i].xyp[j].y2);
	}
	break;
      case FCACT:
	fprintf(data_file," %3.3d_FTYP = %d  %3.3d_FORD = %d  \
%3.3d_REJL = %5.2lf  %3.3d_REJU = %5.2lf  %3.3d_NXYP = %d\n",i+1,act[i].i[0],i+1,
		act[i].i[1],i+1,act[i].d[4],i+1,act[i].d[5],i+1,act[i].nxyp);
	for (j=0; j<act[i].nxyp; j++) {
	  fprintf(data_file,
		  "  %3.3d_XY_%3.3d = %d  %14.8lf  %14.8lf  %14.8lf  %14.8lf\n",
		  i+1,j+1,act[i].xyp[j].i,act[i].xyp[j].x1,act[i].xyp[j].x2,
		  act[i].xyp[j].y1,act[i].xyp[j].y2);
	}
	break;
      case IOACT:
	fprintf(data_file," %3.3d_SPEC = %d  %3.3d_ORDR = %d  %3.3d_NXYP = %d\n",
		i+1,act[i].i[0]+1,i+1,act[i].i[1]+1,i+1,act[i].nxyp);
	for (j=0; j<act[i].nxyp; j++) {
	  fprintf(data_file,"  %3.3d_XY_%3.3d = %14.8lf  %14.8lf\n",i+1,j+1,
		  act[i].xyp[j].x1,act[i].xyp[j].y1);
	}
	break;
      case ICACT:
	fprintf(data_file," %3.3d_NXYP = %d\n",i+1,act[i].nxyp);
	for (j=0; j<act[i].nxyp; j++) {
	  fprintf(data_file,"  %3.3d_XY_%3.3d = %14.8lf  %14.8lf\n",i+1,j+1,
		  act[i].xyp[j].x1,act[i].xyp[j].y1);
	}
	break;
      case SOACT:
	fprintf(data_file," %3.3d_SPEC = %d  %3.3d_ORDR = %d  %3.3d_SCAL =\
 %15.13lf\n",i+1,act[i].i[0]+1,i+1,act[i].i[1]+1,i+1,act[i].d[0]);
	break;
      case ARACT:
	fprintf(data_file," %3.3d_SCLP = %14.8lf  %3.3d_SERR = %14.8lf\n",i+1,
		act[i].d[0],i+1,act[i].d[1]);
	break;
      }
    }
  }

  /* Close file */
  fclose(data_file);

  return 1;

}
