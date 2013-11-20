/******************************************************************************
AST_JD2DATE: Convert Julian date to calendar date

double	j			# Julian day
int	year			# Year
int	month			# Month (1-12)
int	day			# Day of month
double	t			# Time for date (mean solar day)

******************************************************************************/

#include "astron.h"

int ast_jd2date(double j, int *year, int *month, int *day, double *t) {

  int ja=0,jb=0,jc=0,jd=0,je=0;

  ja=AST_NINT(j);
  *t=24.0*(j+0.5-(double)ja);

  if (ja>=2299161) {
    jb=(int)(((double)(ja-1867216)-0.25)/36524.25);
    ja=ja+1+jb-(int)(jb/4);
  }

  jb=ja+1524;
  jc=(int)(6680.0+((double)(jb-2439870)-122.1)/C_JYEAR);
  jd=365*jc+(int)(jc/4);
  je=(int)((double)(jb-jd)/30.6001);
  *day=jb-jd-(int)(30.6001*je);
  *month=je-1;
  if (*month>12) *month-=12;
  *year=jc-4715;
  if (*month>2) (*year)--;
  if (*year<0) (*year)--;

  return 1;

}
