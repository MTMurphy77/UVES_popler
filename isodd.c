/****************************************************************************
Function to decide whether an integer is odd
****************************************************************************/

int isodd(int number) {

  if ((float)(number/2)-((float)number)/2.0<-0.49) return 1;
  return 0;

}
