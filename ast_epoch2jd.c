/******************************************************************************
AST_EPOCH2JD: Convert epoch to Julian Day
******************************************************************************/

#include "astron.h"

double ast_epoch2jd(double epoch) {

  double jd=0.0;

  jd=C_JD2000+(epoch-C_J2000)*C_JYEAR;

  return jd;

}
