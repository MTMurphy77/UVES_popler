/****************************************************************************
* Calculate heliocentric correction for a spectrum and compare with header
* values
****************************************************************************/

#include <math.h>
#include "UVES_popler.h"
#include "astron.h"
#include "error.h"

int UVES_vhelio(spectrum *spec) {

  double   ut=0.0,ra=0.0,dec=0.0;
  double   vrot=0.0,vorb=0.0,vbary=0.0; 
  int      year=0,month=0,day=0;

  /* Determine the observation date and universal time from the Julian Date
     in the header */
  if (!ast_jd2date(spec->jd,&year,&month,&day,&ut))
    errormsg("UVES_vhelio(): Unknown error returned from ast_jd2date()");

  /* Check consistency with date read from header */
  if (year!=spec->year || month!=spec->month || day!=spec->day)
    warnmsg("UVES_vhelio(): Observation date in header (=%4.4d-%2.2d-%2.2d)\n\
\tdoesn't match date calculated from MJD in header (=%4.4d-%2.2d-%2.2d)\n\
\tin file %s.\n\tUsing the latter",
	    spec->year,spec->month,spec->day,year,month,day,spec->file,
	    spec->ut);
  else if (fabs(ut-spec->ut)>UTTOL/60.0) {
    if (spec->ut<UTTOL/60.0 && fabs(ut+spec->ut-24.0)<UTTOL/60.0)
      warnmsg("UVES_vhelio(): UT in header (=%6.3lf hrs) probably\n\
\toccurs next UT-day after recorded obs. date (=%4.4d-%2.2d-%2.2d)\n\
\tin file %s\n\tsince UT calculated from MJD in header is %6.3lf hrs.\n\
\tAssuming MJD gives the most accurate date+time as usual",spec->ut,spec->year,
	      spec->month,spec->day,spec->file,ut);
    else if (spec->ut>24.0-UTTOL/60.0 && fabs(ut+spec->ut-24.0)<UTTOL/60.0)
      warnmsg("UVES_vhelio(): UT in header (=%6.3lf hrs) probably\n\
\toccurs previous UT-day before recorded obs. date (=%4.4d-%2.2d-%2.2d)\n \
\tin file %s\n\tsince UT calculated from MJD in header is %6.3lf hrs.\n	\
\tAssuming MJD gives the most accurate date+time as usual",spec->ut,spec->year,
	      spec->month,spec->day,spec->file,ut);
    else
      warnmsg("UVES_vhelio(): UT in header (=%6.3lf hrs) and UT\n\
\tcalculated from MJD in header (=%6.3lf hrs) differ by more than\n\
\tUTTOL=%5.2lf mins in file\n\t%s.\n\
\tAssuming MJD gives the most accurate date+time",spec->ut,ut,UTTOL,spec->file);
  }
  /* Assume the MJD has the most accurate start-of-exspoure date+time information */
  spec->year=year; spec->month=month; spec->day=day; spec->ut=ut;
  spec->epoch=ast_date2epoch(spec->year,spec->month,spec->day,spec->ut);

  /* Modify Julian day of observation to be approximately in the middle of
     the exposure */
  spec->epoch+=0.5*spec->etime/C_SECYEAR;

  /* Precess the object coordinates (J2000 in most cases) to the epoch in
     which the observations were carried out */
  if (!ast_precess(spec->ra,spec->dec,spec->equ,spec->epoch,&ra,&dec))
    errormsg("UVES_vhelio(): Unknown error returned from ast_precess()");
  /* Calculate the Earth-moon barycentric velocity relative to Sun (annual) */
  vorb=ast_vorbit(ra,dec,spec->epoch);
  /* Calculate the Earth's velocity around Earth-Moon barycentre (lunar) */
  vbary=ast_vbary(ra,dec,spec->epoch);
  /* Calculate the observer's motion relative to the Earth's centre (diurnal)*/
  vrot=ast_vrotate(ra,dec,spec->epoch,spec->lat,spec->lon,spec->alt);
  /* Calculate heliocentric velocity */
  spec->vhel=vorb+vbary+vrot;

  return 1;

}
