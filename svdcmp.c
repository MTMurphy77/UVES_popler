/*************************************************************************** 
SVDCMP: This is the SVDCMP algorithm taken from Numercial Recipes. Here's
what NR has to say about it

Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.V^T. The matrix U replaces a on output. The diagonal
matrix of singular values W is output as a vector w[1..n]. The matrix V (not
the transpose V^T) is output as v[1..n][1..n].

I have modified the routine in the following ways:

0) Double precision and now an integer function
1) Added own error handling
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "fit.h"
#include "memory.h"
#include "error.h"

int svdcmp(double **a, int m, int n, double *w, double **v) {

  int    flag=0,i=0,its=0,j=0,jj=0,k=0,l=0,nm=0;
  double anorm=0.0,c=0.0,f=0.0,g=0.0,h=0.0,s=0.0,scale=0.0,x=0.0,y=0.0,z=0.0;
  double *rv1=NULL;

  /* Allocate memory for relevant arrays */
  if ((rv1=darray(n))==NULL) {
    nferrormsg("svdcmp(): Cannot allocate memory to rv1\n\tarray of size %d",
	       n); return 0;
  }

  for (i=0; i<n; i++) {
    l=i+1; rv1[i]=scale*g; g=s=scale=0.0;
    if (i<m) {
      for (k=i; k<m; k++) scale+=fabs(a[k][i]);
      if (scale) {
	for (k=i; k<m; k++) {
	  a[k][i]/=scale; s+=a[k][i]*a[k][i];
	}
	f=a[i][i]; g=-SIGN(sqrt(s),f); h=f*g-s; a[i][i]=f-g;
	for (j=l; j<n; j++) {
	  for (s=0.0,k=i; k<m; k++) s+=a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i; k<m; k++) a[k][j]+=f*a[k][i];
	}
	for (k=i; k<m; k++) a[k][i]*=scale;
      }
    }
    w[i]=scale*g; g=s=scale=0.0;
    if (i<m && i!=n-1) {
      for (k=l; k<n; k++) scale+=fabs(a[i][k]);
      if (scale) {
	for (k=l; k<n; k++) {
	  a[i][k]/=scale; s+=a[i][k]*a[i][k];
	}
	f=a[i][l]; g=-SIGN(sqrt(s),f); h=f*g-s; a[i][l]=f-g;
	for (k=l; k<n; k++) rv1[k]=a[i][k]/h;
	for (j=l; j<m; j++) {
	  for (s=0.0,k=l; k<n; k++) s+=a[j][k]*a[i][k];
	  for (k=l; k<n; k++) a[j][k]+=s*rv1[k];
	}
	for (k=l; k<n; k++) a[i][k]*=scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1; i>=0; i--) {
    if (i<n-1) {
      if (g) {
	for (j=l; j<n; j++) v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l; j<n; j++) {
	  for (s=0.0,k=l; k<n; k++) s+=a[i][k]*v[k][j];
	  for (k=l; k<n; k++) v[k][j]+=s*v[k][i];
	}
      }
      for (j=l; j<n; j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0; g=rv1[i]; l=i;
  }
  for (i=MIN(m,n)-1; i>=0; i--) {
    l=i+1; g=w[i];
    for (j=l; j<n; j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l; j<n; j++) {
	for (s=0.0,k=l; k<m; k++) s+=a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i; k<m; k++) a[k][j]+=f*a[k][i];
      }
      for (j=i; j<m; j++) a[j][i]*=g;
    }
    else for (j=i; j<m; j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n-1; k>=0; k--) {
    for (its=0; its<SVDNMAX; its++) {
      flag=1;
      for (l=k; l>=0; l--) {
	nm=l-1;
	if (fabs(rv1[l])+anorm==anorm) {
	  flag=0; break;
	}
	if (fabs(w[nm])+anorm==anorm) break;
      }
      if (flag) {
	c=0.0; s=1.0;
	for (i=l; i<=k; i++) {
	  f=s*rv1[i]; rv1[i]=c*rv1[i];
	  if (fabs(f)+anorm==anorm) break;
	  g=w[i]; h=pythag(f,g); w[i]=h;
	  h=1.0/h; c=g*h; s=-f*h;
	  for (j=0; j<m; j++) {
  	    y=a[j][nm]; z=a[j][i]; a[j][nm]=y*c+z*s; a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l==k) {
	if (z<0.0) {
	  w[k]=-z;
	  for (j=0; j<n; j++) v[j][k]=-v[j][k];
	}
	break;
      }
      if (its==29) {
	nferrormsg("svdcmp(): No convergence in %d iterations",SVDNMAX); return 0;
      }
      x=w[l]; nm=k-1; y=w[nm]; g=rv1[nm]; h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y); g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l; j<=nm; j++) {
	i=j+1; g=rv1[i]; y=w[i];
	h=s*g; g=c*g; z=pythag(f,h);
	rv1[j]=z; c=f/z; s=h/z; f=x*c+g*s; g=g*c-x*s; h=y*s; y*=c;
	for (jj=0; jj<n; jj++) {
	  x=v[jj][j]; z=v[jj][i]; v[jj][j]=x*c+z*s; v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h); w[j]=z;
	if (z) {
	  z=1.0/z; c=f*z; s=h*z;
	}
	f=c*g+s*y; x=c*y-s*g;
	for (jj=0; jj<m; jj++) {
	  y=a[jj][j]; z=a[jj][i]; a[jj][j]=y*c+z*s; a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0; rv1[k]=f; w[k]=x;
    }
  }

  /* Clean up */
  free(rv1);
  
  return 1;

}
