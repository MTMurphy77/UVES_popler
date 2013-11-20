/******************************************************************************
AST_DATE2JD: Converts a date (day, month, year) into the Julian day.

int     year                    # Year (The full year AD, e.g. 1812, 2003)
int     month                   # Month (1-12)
int     day                     # Day of month
double  t                       # Time for date (mean solar day)

******************************************************************************/

#include "astron.h"
#include "error.h"

double ast_date2jd(int year, int month, int day, double t) {

  double jd=0.0;
  int    y=0,m=0,d=0;

  y=year;
  if (month<1 || month>12)
    errormsg("ast_date2jd(): Invalid month passed =%d",month);
  else if (month>2) m=month+1;
  else {
    m=month+13; y--;
  }

  jd=(double)((int)(C_JYEAR*y)+(int)(30.6001*m)+day+1720995);
  if (day+31*(m+12*y) >= 588829) {
    d=(int)(y/100); m=(int)(y/400); jd+=(double)(2-d+m);
  }
  jd+=(double)((int)(t*360000.0+0.5))/360000.0/24.0-0.5;

  return jd;

}
