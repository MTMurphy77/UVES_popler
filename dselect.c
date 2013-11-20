/*************************************************************************** 
DSELECT: This is the select algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Returns the kth smallest value in the array arr[1..n]. The input array
will be rearranged to have this value in location arr[k], with all smaller
elements moved to arr[1..k-1] (in arbitrary order) and all larger elements
in arr[k+1..n] (also in arbitrary order)."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Now returns (k+1)th lowest value. That is, if you were to order the
   input array, it will return the value with array index k.
****************************************************************************/

#define DSWAP(a,b) dtmp=(a);(a)=(b);(b)=dtmp;

double dselect(unsigned long k, unsigned long n, double *arr) {

  unsigned long i=0,ir=0,j=0,l=0,mid=0;
  double        a=0.0,dtmp;

  l=0;
  ir=n-1;
  for (;;) {
    if (ir<=l+1) {
      if (ir==l+1 && arr[ir]<arr[l]) { DSWAP(arr[l],arr[ir]); }
      return arr[k];
    }
    else {
      mid=(l+ir) >> 1; i=l+1;
      DSWAP(arr[mid],arr[i]);
      if (arr[i]>arr[ir]) { DSWAP(arr[i],arr[ir]); }
      if (arr[l]>arr[ir]) { DSWAP(arr[l],arr[ir]); }
      if (arr[i]>arr[l]) { DSWAP(arr[i],arr[l]); }
      j=ir; a=arr[l];
      for (;;) {
	do i++; while (arr[i]<a); do j--; while (arr[j]>a);
	if (j<i) break; DSWAP(arr[i],arr[j]);
      }
      arr[l]=arr[j]; arr[j]=a; if (j>=k) ir=j-1; if (j<=k) l=i;
    }
  }
}
