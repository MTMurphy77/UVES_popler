/******************************************************************************
AST_DATE2EPOCH: Convert Gregorian date and solar mean time to a Julian
epoch. A Julian epoch has 365.25 days per year and 24 hours per day.

int	year			# Year
int	month			# Month (1-12)
int	day			# Day of month
double	ut			# Universal time for date (mean solar day)
double	epoch			# Julian epoch

******************************************************************************/

#include "astron.h"

double ast_date2epoch(int year, int month, int day, double ut) {

  double jd=0.0,epoch=0.0;

  jd=ast_date2jd(year,month,day,ut);
  epoch=C_J2000+(jd-C_JD2000)/C_JYEAR;

  return epoch;

}
