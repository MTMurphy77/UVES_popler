/******************************************************************************
AST_JD2EPOCH: Convert Julian Day to epoch
******************************************************************************/

#include "astron.h"

double ast_jd2epoch(double jd) {

  double epoch=0.0;

  epoch=C_J2000+(jd-C_JD2000)/C_JYEAR;

  return epoch;

}
