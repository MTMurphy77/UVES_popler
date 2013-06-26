#!/bin/tcsh

set SMIN = '0.02'
set SMAX = '1.5'
set SSTEP = '0.01'
@ NSTEP = `echo "$SMIN $SMAX $SSTEP" | awk '{printf "%d",1+($2-$1)/$3}'`

/bin/rm -f GFWHM_vs_seeing.dat ; touch GFWHM_vs_seeing.dat

@ i = 0
while ( $i < $NSTEP )
  set S = `echo "$SMIN $SSTEP $i" | awk '{print $1+$2*$3}'`
  awk '{x=0.0; if ($1 > -0.5 && $1 < 0.5) x=exp(-0.5*$1*$1/'$S'/'$S'); printf "%9.6lf %9.6lf 0.100\n",$1,x}' GFWHM_vs_seeing_in.dat > in
  /bin/rm -f funfit.out.1
  set SGUESS = `echo "$S" | awk '{x=0.5; if ($1 < x) x=$1; print x}'`
  printf "n\nin\nn\ng\n0.0\n1.0\n$SGUESS\n0.0\nn\nn\n/NULL\nn\nq\n\n" | funfit > /dev/null
  set GFWHM = `grep "C = " funfit.out.1 | awk '{printf "%8.6lf",$3*2.354820045*1.001*1.006610662}'`
  echo "$S $GFWHM" | awk '{printf "%7.5lf %8.6lf\n",$1*1.001*2.354820045,$2}' >> GFWHM_vs_seeing.dat
  @ i = $i + 1
end

/bin/rm -f in funfit.out.1
