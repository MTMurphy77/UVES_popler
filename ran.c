/*************************************************************************** 
RAN: This is the ran2 algorithm taken from Numercial Recipes. The period is
(apparently) practically inexhaustable, about 2.3x10^(18). Serial
correlations are avoided using two methods: combining two generators and by
using a shuffle of the output numbers. NR offer $1000 to anyone who can
show this routine not to give random numbers. Here's what NR has to say
about it:

Long period (>2x10^(18)) random number generator of L'Ecuyer with
Bays-Durham shuffle and added safeguards. Returns a uniform random 
deviate between 0.0 and 1.0 (exlcusive of the endpoint values). Call 
with idum a negative integer to initialize; thereafter, do no alter 
idum between successive deviates in a sequence. RNMX should approximate 
the largest floating value that is less than 1.
****************************************************************************/

#define IM1  2147483563
#define IM2  2147483399
#define AM   (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1  40014
#define IA2  40692
#define IQ1  53668
#define IQ2  52774
#define IR1  12211
#define IR2  3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS  1.2e-7
#define RNMX (1.0-EPS)

double ran(long *idum) {

  int         j;
  long        k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (*idum<=0) {
    if (-(*idum)<1) *idum=1;
    else *idum=-(*idum);
    idum2=(*idum);
    for (j=NTAB+7; j>=0; j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum<0) *idum+=IM1;
      if (j<NTAB) iv[j]=*idum;
    }
    iy=iv[0];
  }

  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum<0) *idum+=IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2<0) idum2+=IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j]=*idum;
  if (iy<1) iy+=IMM1;
  if ((temp=AM*iy)>RNMX) return RNMX;
  else return temp;

}
