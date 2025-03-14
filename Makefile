# Makefile for UVES_popler

# Linux
# NOTE: Change compilation command for "UVES_popler" below to use
# $(FC) for compilation and linking as well.
#SHELL = tcsh
#CC = gcc
#FC = gfortran
#CFLAGS = -L${CFITSIO_DIR}/include -L${PGPLOT_DIR} -O2 -Wall -I./
#LIBS = -L${CFITSIO_DIR}/lib -L${PGPLOT_DIR} -lm -lcpgplot -lpgplot -lX11 -lcfitsio

# Mac OS X - assumes PGPLOT and CFITSIO were installed through Homebrew.
SHELL = tcsh
CC = gcc -Os
FC = gfortran -Os
CFLAGS = -O2 -Wall -I./ -I${PGPLOT_DIR}
LIBS = -ldl -lm -L/opt/homebrew/lib -lpng -lcurl -L${PGPLOT_DIR} -lcpgplot -lpgplot -lX11 -lcfitsio

TARGET = ${HOME}/bin/UVES_popler
DEVTARGET = ${HOME}/bin/UVES_popler.dev

OBJECTS = ast_coord.o ast_date2epoch.o ast_date2jd.o ast_epoch2jd.o ast_jd2date.o ast_jd2epoch.o ast_mst.o ast_precess.o ast_rotmatrix.o ast_vbary.o ast_vorbit.o ast_vrotate.o barray.o brent.o carray.o chebyshev_eval.o chixy.o chixy_fixa.o cmatrix.o covsrt.o dpolint.o dselect.o edlen_a2v.o erffn.o errormsg.o EW.o darray.o djmax.o djmin.o dmatrix.o farray.o faskropen.o faskwopen.o fcompl.o fitexy.o fjmax.o fjmin.o gammcf.o gammln.o gammp.o gammq.o gammser.o gaussj.o get_input.o getscbc.o iarray.o idxdmax.o idxdval.o idxfval.o idxival.o ijmax.o imatrix.o isodd.o isdir.o legendre_eval.o linfit.o median.o medianrun.o mnbrak.o mrqcof.o mrqfit_erffn.o mrqfit_gauss.o mrqfit_multierffn.o mrqmin.o nferrormsg.o pg_button.o pg_get_wins.o pg_open.o pg_win_rename.o poly_eval.o pythag.o qsort_darray.o qsort_dbleint.o qsort_dbletwointarray.o ran.o sigclip.o spline.o splint.o stats.o strisnum.o strlower.o svbksb.o svdcmp.o svdfit.o svdfit_chebyshev.o svdfit_legendre.o svdfit_poly.o svdvar.o UVES_atmask.o UVES_boxcar.o UVES_chunk_cont.o UVES_combine_cont.o UVES_combine_nocont.o UVES_combine_region.o UVES_combine_spec.o UVES_combsynthThAr.o UVES_confit.o UVES_cspec_cont.o UVES_cspec_stats.o UVES_hex2dec.o UVES_hirx_blzfit.o UVES_init_cspec.o UVES_memspec.o UVES_merge_thar.o UVES_model_resol.o UVES_order_cont.o UVES_order_rejsigedge.o UVES_order_sigclip.o UVES_order_stats.o UVES_params_init.o UVES_params_set.o UVES_past_actions.o UVES_pgenv_init.o UVES_pixscal.o UVES_plot_cspec.o UVES_replace_envinstr.o UVES_plot_replay.o UVES_popler.o UVES_r1Dspec.o UVES_r2Dspec.o UVES_r2Dspec_ESOmer.o UVES_r2Dspec_espresso.o UVES_r2Dspec_harps.o UVES_r2Dspec_hirx.o UVES_r2Dspec_iraf.o UVES_r2Dspec_iresi.o UVES_r2Dspec_irls.o UVES_r2Dspec_KODIAQ.o UVES_r2Dspec_mage.o UVES_r2Dspec_makee.o UVES_r2Dspec_pypeit.o UVES_ratmask.o UVES_redispers.o UVES_replay_control.o UVES_rescale_region.o UVES_revwpol.o UVES_rinputfile.o UVES_rFITSlist.o UVES_rMacmap.o UVES_rscale.o UVES_rUPLfile.o UVES_rvshift.o UVES_scale.o UVES_select_subspec.o UVES_set_wavelen_scale.o UVES_skysub.o UVES_synthThAr.o UVES_thar_sigarray.o UVES_undo_lastact.o UVES_vhelio.o UVES_wDATfile.o UVES_wFITSfile.o UVES_wrawFITS.o UVES_wUPLfile.o UVES_wpol.o UVES_tmpfit.o warnmsg.o zbrent.o

UVES_popler: $(OBJECTS)
#	$(CC) -o $(TARGET) $(OBJECTS) $(LIBS)
	$(FC) -o $(TARGET) $(OBJECTS) $(LIBS)

dev: $(OBJECTS)
#	$(CC) -o $(DEVTARGET) $(OBJECTS) $(LIBS)
	$(FC) -o $(DEVTARGET) $(OBJECTS) $(LIBS)

depend:
	makedepend -f Makefile -Y -- $(CFLAGS) -- -s "# Dependencies" \
	$(OBJECTS:.o=.c) $(FORT-OBJECTS:.o=.f) >& /dev/null

clean: 
	/bin/rm -f *~ *.o

# Dependencies

ast_date2epoch.o: astron.h const.h
ast_date2jd.o: astron.h const.h error.h
ast_epoch2jd.o: astron.h const.h
ast_jd2date.o: astron.h const.h
ast_jd2epoch.o: astron.h const.h
ast_mst.o: astron.h const.h
ast_precess.o: astron.h const.h memory.h error.h
ast_rotmatrix.o: astron.h const.h
ast_vbary.o: astron.h const.h error.h
ast_vorbit.o: astron.h const.h error.h
ast_vrotate.o: astron.h const.h
barray.o: error.h
brent.o: fit.h const.h error.h
carray.o: error.h
chebyshev_eval.o: error.h
chixy.o: fit.h
chixy_fixa.o: fit.h
cmatrix.o: memory.h error.h
covsrt.o: fit.h
dpolint.o: memory.h error.h
edlen_a2v.o: error.h
erffn.o: gamm.h
EW.o: fit.h gamm.h memory.h utils.h const.h error.h
darray.o: error.h
djmax.o: error.h
djmin.o: error.h
dmatrix.o: memory.h error.h
farray.o: error.h
faskropen.o: file.h input.h error.h
faskwopen.o: file.h input.h error.h
fcompl.o: charstr.h file.h error.h
fitexy.o: fit.h memory.h stats.h const.h error.h
fjmax.o: error.h
fjmin.o: error.h
gammcf.o: gamm.h error.h
gammln.o: gamm.h
gammp.o: gamm.h error.h
gammq.o: gamm.h error.h
gammser.o: gamm.h error.h
gaussj.o: fit.h memory.h error.h
get_input.o: charstr.h error.h input.h
getscbc.o: charstr.h input.h error.h
iarray.o: error.h
ijmax.o: error.h
imatrix.o: memory.h error.h
legendre_eval.o: fit.h
linfit.o: error.h
median.o: sort.h stats.h memory.h error.h
medianrun.o: stats.h sort.h memory.h error.h
mnbrak.o: fit.h const.h
mrqcof.o: fit.h memory.h error.h
mrqfit_erffn.o: fit.h gamm.h memory.h const.h error.h
mrqfit_gauss.o: fit.h memory.h error.h
mrqfit_multierffn.o: fit.h gamm.h memory.h const.h error.h
mrqmin.o: fit.h memory.h error.h
pg_button.o: pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h utils.h
pg_get_wins.o: pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h utils.h
pg_get_wins.o: error.h
pg_open.o: pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h utils.h error.h
pg_win_rename.o: pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h utils.h
pg_win_rename.o: error.h
pythag.o: fit.h
qsort_dbleint.o: sort.h
qsort_dbletwointarray.o: sort.h
sigclip.o: stats.h memory.h error.h
spline.o: memory.h error.h
splint.o: error.h
stats.o: stats.h error.h
strisnum.o: charstr.h error.h
svbksb.o: memory.h error.h
svdcmp.o: fit.h memory.h error.h
svdfit.o: fit.h memory.h error.h
svdfit_chebyshev.o: error.h
svdvar.o: memory.h error.h
UVES_atmask.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_atmask.o: utils.h memory.h error.h const.h
UVES_boxcar.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_boxcar.o: utils.h stats.h const.h memory.h error.h
UVES_chunk_cont.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_chunk_cont.o: charstr.h utils.h const.h memory.h error.h
UVES_combine_cont.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_combine_cont.o: charstr.h utils.h stats.h memory.h error.h
UVES_combine_nocont.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_combine_nocont.o: charstr.h utils.h stats.h memory.h error.h
UVES_combine_region.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_combine_region.o: charstr.h utils.h stats.h fit.h memory.h error.h
UVES_combine_spec.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_combine_spec.o: charstr.h utils.h stats.h memory.h error.h
UVES_combsynthThAr.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_combsynthThAr.o: charstr.h utils.h gamm.h const.h error.h
UVES_confit.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_confit.o: utils.h fit.h sort.h memory.h error.h
UVES_cspec_cont.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_cspec_cont.o: charstr.h utils.h stats.h const.h error.h
UVES_cspec_stats.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_cspec_stats.o: charstr.h utils.h memory.h stats.h const.h error.h
UVES_hex2dec.o: error.h
UVES_hirx_blzfit.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_hirx_blzfit.o: charstr.h utils.h fit.h stats.h memory.h error.h
UVES_init_cspec.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_init_cspec.o: charstr.h utils.h const.h memory.h error.h
UVES_memspec.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_memspec.o: utils.h memory.h error.h
UVES_merge_thar.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_merge_thar.o: charstr.h utils.h const.h memory.h error.h
UVES_model_resol.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_model_resol.o: charstr.h utils.h fit.h stats.h memory.h const.h error.h
UVES_order_cont.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_order_cont.o: charstr.h utils.h error.h
UVES_order_rejsigedge.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_order_rejsigedge.o: charstr.h utils.h memory.h stats.h error.h
UVES_order_sigclip.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_order_sigclip.o: charstr.h utils.h stats.h memory.h error.h
UVES_order_stats.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_order_stats.o: charstr.h utils.h stats.h memory.h error.h
UVES_params_init.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_params_init.o: charstr.h utils.h
UVES_params_set.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_params_set.o: charstr.h utils.h
UVES_past_actions.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_past_actions.o: charstr.h utils.h fit.h memory.h error.h
UVES_pgenv_init.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_pgenv_init.o: charstr.h utils.h
UVES_pixscal.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_pixscal.o: utils.h const.h
UVES_plot_cspec.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_plot_cspec.o: charstr.h utils.h fit.h stats.h sort.h gamm.h input.h
UVES_plot_cspec.o: memory.h const.h error.h
UVES_replace_envinstr.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_replace_envinstr.o: charstr.h utils.h memory.h error.h
UVES_plot_replay.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_plot_replay.o: charstr.h utils.h input.h error.h
UVES_popler.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_popler.o: utils.h stats.h input.h error.h
UVES_r1Dspec.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_r1Dspec.o: utils.h const.h file.h memory.h error.h
UVES_r2Dspec.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_r2Dspec.o: utils.h stats.h memory.h error.h const.h
UVES_r2Dspec_ESOmer.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_ESOmer.o: charstr.h utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_espresso.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_espresso.o: charstr.h utils.h stats.h memory.h error.h const.h
UVES_r2Dspec_harps.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_harps.o: charstr.h utils.h stats.h memory.h error.h const.h
UVES_r2Dspec_hirx.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_hirx.o: charstr.h utils.h fit.h memory.h error.h const.h
UVES_r2Dspec_iraf.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_iraf.o: charstr.h utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_iresi.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_iresi.o: charstr.h utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_irls.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_irls.o: charstr.h utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_KODIAQ.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_KODIAQ.o: charstr.h utils.h memory.h error.h
UVES_r2Dspec_mage.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_mage.o: charstr.h utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_makee.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_makee.o: charstr.h utils.h memory.h const.h error.h
UVES_r2Dspec_pypeit.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_r2Dspec_pypeit.o: charstr.h utils.h fit.h memory.h error.h const.h
UVES_ratmask.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_ratmask.o: utils.h memory.h file.h error.h
UVES_redispers.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_redispers.o: charstr.h utils.h memory.h error.h
UVES_replay_control.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_replay_control.o: charstr.h utils.h input.h stats.h memory.h error.h
UVES_replay_control.o: const.h
UVES_rescale_region.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_rescale_region.o: charstr.h utils.h stats.h fit.h memory.h error.h
UVES_revwpol.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_revwpol.o: utils.h fit.h memory.h error.h
UVES_rinputfile.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_rinputfile.o: charstr.h utils.h file.h error.h
UVES_rFITSlist.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_rFITSlist.o: charstr.h utils.h file.h error.h
UVES_rMacmap.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_rMacmap.o: utils.h memory.h file.h error.h
UVES_rscale.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_rscale.o: utils.h file.h error.h
UVES_rUPLfile.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_rUPLfile.o: charstr.h utils.h file.h error.h
UVES_rvshift.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_rvshift.o: utils.h file.h error.h
UVES_scale.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_scale.o: utils.h error.h
UVES_select_subspec.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_select_subspec.o: charstr.h utils.h memory.h error.h
UVES_set_wavelen_scale.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_set_wavelen_scale.o: charstr.h utils.h const.h memory.h error.h
UVES_skysub.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_skysub.o: utils.h memory.h error.h
UVES_synthThAr.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_synthThAr.o: charstr.h utils.h gamm.h const.h error.h
UVES_thar_sigarray.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_thar_sigarray.o: charstr.h utils.h memory.h error.h
UVES_undo_lastact.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_undo_lastact.o: charstr.h utils.h error.h
UVES_vhelio.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_vhelio.o: utils.h astron.h const.h error.h
UVES_wDATfile.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_wDATfile.o: charstr.h utils.h file.h error.h
UVES_wFITSfile.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_wFITSfile.o: charstr.h utils.h memory.h error.h
UVES_wrawFITS.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_wrawFITS.o: charstr.h utils.h memory.h error.h
UVES_wUPLfile.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h
UVES_wUPLfile.o: charstr.h utils.h file.h error.h
UVES_wpol.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_wpol.o: utils.h fit.h stats.h memory.h const.h error.h
UVES_tmpfit.o: UVES_popler.h pg_plot.h /usr/local/pgplot/cpgplot.h charstr.h
UVES_tmpfit.o: utils.h fit.h memory.h const.h error.h
zbrent.o: fit.h error.h
