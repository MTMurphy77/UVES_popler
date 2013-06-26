/****************************************************************************
* Break a given piece of spectrum up into chunks of specified size (in
* velocity), fit continua to each chunk and then join the continuum chunks
****************************************************************************/

#include <stdlib.h>
#include "UVES_popler.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int UVES_chunk_cont(double *wl, double *fl, double *er, double *co, int *st,
		    int n, double vc, double rsigl, double rsigu, double pctl,
		    int typ, int ord) {

  double   weight=0.0,dweight=0.0;
  double   ddum=0.0;
  double   *fit=NULL;
  int      nfitsec=0,fit_nmax=0,fit_n=0,fit_nl=0,fit_nr=0;
  int      sidx=0;
  int      *fit_nsec=NULL,**fit_idxsec=NULL;
  int      i=0,j=0,k=0,idum=0;

  /* Decide how many fitting regions are contained between start and end */
  /* The following assumes that the centre of the first section is defined
     with respect to the start of the data array, all other centres being
     defined with respect to the previous centre. The ends of each section
     are defined with respect to the next centre. The buffer starts are
     defined as the centres of the previouss section and the buffer ends
     are defined as the centres of the following sections */
  nfitsec=1; sidx=0;
  if ((fit_nsec=iarray(nfitsec))==NULL)
    errormsg("UVES_chunk_cont(): Cannot allocate memory to fit_nsec\n\
\tarray of size %d",nfitsec);
  idum=((idum=idxdval(wl,n,wl[sidx]*(1.0+0.5*vc/C_C_K))-1)>=0) ? idum : n-1;
  fit_nsec[nfitsec-1]=idum-sidx+1; sidx=idum;
  if (idum<n-1-ord-1) {
    ddum=wl[idum];
    while ((idum=idxdval(wl,n,ddum*(1.0+vc/C_C_K))-1)>=0 && idum<=n-1-ord-1) {
      ddum=wl[idum]; nfitsec++;
      if (!(fit_nsec=(int *)realloc(fit_nsec,(size_t)(nfitsec*sizeof(int)))))
	errormsg("UVES_chunk_cont(): Cannot reallocate memory to fit_nsec\n\
\tarray of size %d",nfitsec); 
      fit_nsec[nfitsec-1]=idum-sidx+1; sidx=idum;
    }
  }
  
  /* Allocate memory for fit index matrix */
  if ((fit_idxsec=imatrix(nfitsec,5))==NULL)
    errormsg("UVES_chunk_cont(): Cannot allocate memory to fit_idxsec\n\
\tarray of size %dx%d",nfitsec,5);
  
  /* Fill fit_idxsec array with indices of fitting sections and buffers */
  for (i=0,fit_nmax=0; i<nfitsec; i++) {
    if (!i) {
      fit_idxsec[i][0]=fit_idxsec[i][1]=0;
      if (nfitsec==1) {
	fit_idxsec[i][2]=fit_idxsec[i][3]=fit_idxsec[i][4]=n-1;
	fit_nmax=fit_idxsec[i][4]-fit_idxsec[i][0]+1;
      }
      else fit_idxsec[i][2]=fit_idxsec[i][0]+fit_nsec[0]-1;
    }
    else if (i==nfitsec-1) {
      fit_idxsec[i][0]=fit_idxsec[i-1][2];
      fit_idxsec[i][3]=fit_idxsec[i][4]=n-1;
      fit_idxsec[i][2]=fit_idxsec[i-1][4]=fit_idxsec[i][0]+fit_nsec[i]-1;
      fit_idxsec[i][1]=fit_idxsec[i-1][3]=fit_idxsec[i][2]-fit_nsec[i]/2;
      fit_nmax=(MAX(fit_nmax,(fit_idxsec[i-1][4]-fit_idxsec[i-1][0]+1)));
      fit_nmax=(MAX(fit_nmax,(fit_idxsec[i][4]-fit_idxsec[i][0]+1)));
    }
    else {
      fit_idxsec[i][0]=fit_idxsec[i-1][2];
      fit_idxsec[i][2]=fit_idxsec[i-1][4]=fit_idxsec[i][0]+fit_nsec[i]-1;
      fit_idxsec[i][1]=fit_idxsec[i-1][3]=fit_idxsec[i][2]-fit_nsec[i]/2;
      fit_nmax=(MAX(fit_nmax,(fit_idxsec[i-1][4]-fit_idxsec[i-1][0]+1)));
    }
  }
  
  /* Allocate memory for fitting array */
  if ((fit=darray(fit_nmax))==NULL)
      errormsg("UVES_chunk_cont(): Cannot allocate memory to fit\n\
\tarray of size %d",fit_nmax);
  
  /* Fit first section+buffer */
  fit_n=fit_idxsec[0][4]-fit_idxsec[0][0]+1;
  idum=UVES_confit(&(fl[fit_idxsec[0][0]]),&(er[fit_idxsec[0][0]]),
		   &(st[fit_idxsec[0][0]]),fit_n,typ,ord,rsigl,rsigu,pctl,0,fit);
  if (!idum) {
    nferrormsg("UVES_chunk_cont(): Unable to fit continuum to\n\
\tcombined spectrum between pixels %d and %d,\n\
\twavelengths %lf and %lf, in UVES_confit()",fit_idxsec[0][0],fit_idxsec[0][4],
	     wl[fit_idxsec[0][0]],wl[fit_idxsec[0][4]]); return 0;
  } else if (idum==1)
    for (i=0,j=fit_idxsec[0][0]; i<fit_n; i++,j++) co[j]=fit[i];

  /* Now each consecutive fitting section+buffer and merge them with
     existing continuuim */
  for (i=1; i<nfitsec; i++) {
    fit_n=fit_idxsec[i][4]-fit_idxsec[i][0]+1;
    idum=UVES_confit(&(fl[fit_idxsec[i][0]]),&(er[fit_idxsec[i][0]]),
		     &(st[fit_idxsec[i][0]]),fit_n,typ,ord,rsigl,rsigu,pctl,0,fit);
    if (!idum) {
      nferrormsg("UVES_chunk_cont(): Unable to fit continuum to\n\
\tcombined spectrum between pixels %d and %d,\n\
\twavelengths %lf and %lf, in UVES_confit()",fit_idxsec[i][0],fit_idxsec[i][4],
	       wl[fit_idxsec[i][0]],wl[fit_idxsec[i][4]]); return 0;
    } else if (idum==1) {
      fit_nl=fit_idxsec[i][1]-fit_idxsec[i][0]+1;
      fit_nr=fit_idxsec[i][2]-fit_idxsec[i][1]+1;
      if (fit_nl<=1 || fit_nr<=1) {
	nferrormsg("UVES_chunk_cont(): Half-length of continuum fitting\n\
\tsection is too small. Left and right half-lengths are %d and %d pixels",
		   fit_nl,fit_nr); return 0;
      }
      dweight=1.0/(double)(fit_nl-1);
      for (j=1,k=fit_idxsec[i][0]+1,weight=dweight; j<fit_nl; j++,k++) {
	co[k]+=weight*fit[j]; co[k]/=1.0+weight; weight+=dweight;
      }
      dweight=1.0/(double)(fit_nr-1);
      for (j=fit_nl,k=fit_idxsec[i][1]+1,weight=1.0-dweight; j<fit_nl+fit_nr-2;
	   j++,k++) {
	co[k]=weight*co[k]+fit[j]; co[k]/=1.0+weight;
	weight-=dweight;
      }
      for (j=fit_nl+fit_nr-2,k=fit_idxsec[i][2]; k<=fit_idxsec[i][4]; j++,k++)
	co[k]=fit[j];
    }
  }

  /* Clean up */
  free(fit); free(fit_nsec); free(*fit_idxsec); free(fit_idxsec);
  
  return 1;

}
