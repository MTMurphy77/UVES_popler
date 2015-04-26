# Makefile for UVES_popler

# Linux
# NOTE: Change compilation command for "UVES_popler" below to use
# $(FC) for compilation and linking as well.
#SHELL = tcsh
#CC = gcc
#FC = gfortran
#CFLAGS = -O2 -Wall -I./
#LIBS += -lm -lcpgplot -lpgplot -lX11 -lcfitsio

# Mac OS X - assumes PGPLOT and CFITSIO were installed through MacPorts.
SHELL = tcsh
CC = gcc
CFLAGS = -O2 -Wall -I./ -I/opt/local/include
LIBS = -ldl -lm -lX11 -L/opt/local/lib -lcpgplot -lpgplot /opt/local/lib/libcfitsio.a

TARGET = ${HOME}/bin/UVES_popler
DEVTARGET = ${HOME}/bin/UVES_popler.dev

OBJECTS = ast_coord.o ast_date2epoch.o ast_date2jd.o ast_epoch2jd.o ast_jd2date.o ast_jd2epoch.o ast_mst.o ast_precess.o ast_rotmatrix.o ast_vbary.o ast_vorbit.o ast_vrotate.o brent.o carray.o chebyshev_eval.o chixy.o chixy_fixa.o cmatrix.o covsrt.o dpolint.o dselect.o edlen_a2v.o edlen_v2a.o erffn.o errormsg.o EW.o darray.o djmax.o djmin.o dmatrix.o farray.o faskropen.o faskwopen.o fcompl.o fitexy.o fjmax.o fjmin.o gammcf.o gammln.o gammp.o gammq.o gammser.o gaussj.o get_input.o getscbc.o iarray.o idxdmax.o idxdval.o idxfval.o idxival.o ijmax.o imatrix.o isodd.o isdir.o legendre_eval.o linfit.o median.o medianrun.o mnbrak.o mrqcof.o mrqfit_erffn.o mrqfit_gauss.o mrqfit_multierffn.o mrqmin.o nferrormsg.o pg_button.o pg_get_wins.o pg_open.o pg_win_rename.o poly_eval.o pythag.o qsort_darray.o qsort_dbleint.o qsort_dbletwointarray.o ran.o sigclip.o spline.o splint.o stats.o strisnum.o strlower.o svbksb.o svdcmp.o svdfit.o svdfit_chebyshev.o svdfit_legendre.o svdfit_poly.o svdvar.o UVES_atmask.o UVES_boxcar.o UVES_chunk_cont.o UVES_combine_cont.o UVES_combine_nocont.o UVES_combine_region.o UVES_combine_spec.o UVES_combsynthThAr.o UVES_confit.o UVES_cspec_cont.o UVES_hex2dec.o UVES_hirx_blzfit.o UVES_init_cspec.o UVES_memspec.o UVES_merge_thar.o UVES_model_resol.o UVES_order_cont.o UVES_order_rejsigedge.o UVES_order_sigclip.o UVES_order_stats.o UVES_params_init.o UVES_params_set.o UVES_past_actions.o UVES_pgenv_init.o UVES_pixscal.o UVES_plot_cspec.o UVES_plot_replay.o UVES_popler.o UVES_r1Dspec.o UVES_r2Dspec.o UVES_r2Dspec_ESOmer.o UVES_r2Dspec_harps.o UVES_r2Dspec_hirx.o UVES_r2Dspec_iraf.o UVES_r2Dspec_iresi.o UVES_r2Dspec_irls.o UVES_r2Dspec_mage.o UVES_r2Dspec_makee.o UVES_ratmask.o UVES_redispers.o UVES_replay_control.o UVES_rescale_region.o UVES_revwpol.o UVES_rinputfile.o UVES_rFITSlist.o UVES_rMacmap.o UVES_rUPLfile.o UVES_rvshift.o UVES_select_subspec.o UVES_set_wavelen_scale.o UVES_skysub.o UVES_synthThAr.o UVES_thar_sigarray.o UVES_undo_lastact.o UVES_vhelio.o UVES_wDATfile.o UVES_wFITSfile.o UVES_wrawFITS.o UVES_wUPLfile.o UVES_wpol.o UVES_tmpfit.o warnmsg.o zbrent.o

UVES_popler: $(OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS) $(LIBS)
#	$(FC) -o $(TARGET) $(OBJECTS) $(LIBS)

dev: $(OBJECTS)
	$(CC) -o $(DEVTARGET) $(OBJECTS) $(LIBS)
#	$(FC) -o $(DEVTARGET) $(OBJECTS) $(LIBS)

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
brent.o: fit.h const.h error.h
carray.o: error.h
chebyshev_eval.o: error.h
chixy.o: fit.h
chixy_fixa.o: fit.h
cmatrix.o: memory.h error.h
covsrt.o: fit.h
dpolint.o: memory.h error.h
edlen_a2v.o: error.h
edlen_v2a.o: error.h
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
pg_button.o: pg_plot.h /opt/local/include/cpgplot.h
pg_button.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
pg_button.o: /opt/local/include/X11/Xfuncproto.h
pg_button.o: /opt/local/include/X11/Xosdefs.h /opt/local/include/X11/Xutil.h
pg_button.o: /opt/local/include/X11/keysym.h
pg_button.o: /opt/local/include/X11/keysymdef.h
pg_button.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
pg_button.o: /opt/local/include/X11/Xarch.h
pg_button.o: /opt/local/include/X11/extensions/shape.h
pg_button.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h utils.h
pg_get_wins.o: pg_plot.h /opt/local/include/cpgplot.h
pg_get_wins.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
pg_get_wins.o: /opt/local/include/X11/Xfuncproto.h
pg_get_wins.o: /opt/local/include/X11/Xosdefs.h
pg_get_wins.o: /opt/local/include/X11/Xutil.h /opt/local/include/X11/keysym.h
pg_get_wins.o: /opt/local/include/X11/keysymdef.h
pg_get_wins.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
pg_get_wins.o: /opt/local/include/X11/Xarch.h
pg_get_wins.o: /opt/local/include/X11/extensions/shape.h
pg_get_wins.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
pg_get_wins.o: utils.h error.h
pg_open.o: pg_plot.h /opt/local/include/cpgplot.h
pg_open.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
pg_open.o: /opt/local/include/X11/Xfuncproto.h
pg_open.o: /opt/local/include/X11/Xosdefs.h /opt/local/include/X11/Xutil.h
pg_open.o: /opt/local/include/X11/keysym.h /opt/local/include/X11/keysymdef.h
pg_open.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
pg_open.o: /opt/local/include/X11/Xarch.h
pg_open.o: /opt/local/include/X11/extensions/shape.h
pg_open.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h utils.h
pg_open.o: error.h
pg_win_rename.o: pg_plot.h /opt/local/include/cpgplot.h
pg_win_rename.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
pg_win_rename.o: /opt/local/include/X11/Xfuncproto.h
pg_win_rename.o: /opt/local/include/X11/Xosdefs.h
pg_win_rename.o: /opt/local/include/X11/Xutil.h
pg_win_rename.o: /opt/local/include/X11/keysym.h
pg_win_rename.o: /opt/local/include/X11/keysymdef.h
pg_win_rename.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
pg_win_rename.o: /opt/local/include/X11/Xarch.h
pg_win_rename.o: /opt/local/include/X11/extensions/shape.h
pg_win_rename.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
pg_win_rename.o: utils.h error.h
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
UVES_atmask.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_atmask.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_atmask.o: /opt/local/include/X11/Xfuncproto.h
UVES_atmask.o: /opt/local/include/X11/Xosdefs.h
UVES_atmask.o: /opt/local/include/X11/Xutil.h /opt/local/include/X11/keysym.h
UVES_atmask.o: /opt/local/include/X11/keysymdef.h
UVES_atmask.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_atmask.o: /opt/local/include/X11/Xarch.h
UVES_atmask.o: /opt/local/include/X11/extensions/shape.h
UVES_atmask.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_atmask.o: utils.h memory.h error.h const.h
UVES_boxcar.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_boxcar.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_boxcar.o: /opt/local/include/X11/Xfuncproto.h
UVES_boxcar.o: /opt/local/include/X11/Xosdefs.h
UVES_boxcar.o: /opt/local/include/X11/Xutil.h /opt/local/include/X11/keysym.h
UVES_boxcar.o: /opt/local/include/X11/keysymdef.h
UVES_boxcar.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_boxcar.o: /opt/local/include/X11/Xarch.h
UVES_boxcar.o: /opt/local/include/X11/extensions/shape.h
UVES_boxcar.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_boxcar.o: utils.h stats.h const.h memory.h error.h
UVES_chunk_cont.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_chunk_cont.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_chunk_cont.o: /opt/local/include/X11/Xfuncproto.h
UVES_chunk_cont.o: /opt/local/include/X11/Xosdefs.h
UVES_chunk_cont.o: /opt/local/include/X11/Xutil.h
UVES_chunk_cont.o: /opt/local/include/X11/keysym.h
UVES_chunk_cont.o: /opt/local/include/X11/keysymdef.h
UVES_chunk_cont.o: /opt/local/include/X11/Xatom.h
UVES_chunk_cont.o: /opt/local/include/X11/Xos.h
UVES_chunk_cont.o: /opt/local/include/X11/Xarch.h
UVES_chunk_cont.o: /opt/local/include/X11/extensions/shape.h
UVES_chunk_cont.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_chunk_cont.o: utils.h const.h memory.h error.h
UVES_combine_cont.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_combine_cont.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_combine_cont.o: /opt/local/include/X11/Xfuncproto.h
UVES_combine_cont.o: /opt/local/include/X11/Xosdefs.h
UVES_combine_cont.o: /opt/local/include/X11/Xutil.h
UVES_combine_cont.o: /opt/local/include/X11/keysym.h
UVES_combine_cont.o: /opt/local/include/X11/keysymdef.h
UVES_combine_cont.o: /opt/local/include/X11/Xatom.h
UVES_combine_cont.o: /opt/local/include/X11/Xos.h
UVES_combine_cont.o: /opt/local/include/X11/Xarch.h
UVES_combine_cont.o: /opt/local/include/X11/extensions/shape.h
UVES_combine_cont.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_combine_cont.o: utils.h stats.h memory.h error.h
UVES_combine_nocont.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_combine_nocont.o: /opt/local/include/X11/Xlib.h
UVES_combine_nocont.o: /opt/local/include/X11/X.h
UVES_combine_nocont.o: /opt/local/include/X11/Xfuncproto.h
UVES_combine_nocont.o: /opt/local/include/X11/Xosdefs.h
UVES_combine_nocont.o: /opt/local/include/X11/Xutil.h
UVES_combine_nocont.o: /opt/local/include/X11/keysym.h
UVES_combine_nocont.o: /opt/local/include/X11/keysymdef.h
UVES_combine_nocont.o: /opt/local/include/X11/Xatom.h
UVES_combine_nocont.o: /opt/local/include/X11/Xos.h
UVES_combine_nocont.o: /opt/local/include/X11/Xarch.h
UVES_combine_nocont.o: /opt/local/include/X11/extensions/shape.h
UVES_combine_nocont.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_combine_nocont.o: charstr.h utils.h stats.h memory.h error.h
UVES_combine_region.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_combine_region.o: /opt/local/include/X11/Xlib.h
UVES_combine_region.o: /opt/local/include/X11/X.h
UVES_combine_region.o: /opt/local/include/X11/Xfuncproto.h
UVES_combine_region.o: /opt/local/include/X11/Xosdefs.h
UVES_combine_region.o: /opt/local/include/X11/Xutil.h
UVES_combine_region.o: /opt/local/include/X11/keysym.h
UVES_combine_region.o: /opt/local/include/X11/keysymdef.h
UVES_combine_region.o: /opt/local/include/X11/Xatom.h
UVES_combine_region.o: /opt/local/include/X11/Xos.h
UVES_combine_region.o: /opt/local/include/X11/Xarch.h
UVES_combine_region.o: /opt/local/include/X11/extensions/shape.h
UVES_combine_region.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_combine_region.o: charstr.h utils.h stats.h fit.h memory.h error.h
UVES_combine_spec.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_combine_spec.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_combine_spec.o: /opt/local/include/X11/Xfuncproto.h
UVES_combine_spec.o: /opt/local/include/X11/Xosdefs.h
UVES_combine_spec.o: /opt/local/include/X11/Xutil.h
UVES_combine_spec.o: /opt/local/include/X11/keysym.h
UVES_combine_spec.o: /opt/local/include/X11/keysymdef.h
UVES_combine_spec.o: /opt/local/include/X11/Xatom.h
UVES_combine_spec.o: /opt/local/include/X11/Xos.h
UVES_combine_spec.o: /opt/local/include/X11/Xarch.h
UVES_combine_spec.o: /opt/local/include/X11/extensions/shape.h
UVES_combine_spec.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_combine_spec.o: utils.h stats.h memory.h error.h
UVES_combsynthThAr.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_combsynthThAr.o: /opt/local/include/X11/Xlib.h
UVES_combsynthThAr.o: /opt/local/include/X11/X.h
UVES_combsynthThAr.o: /opt/local/include/X11/Xfuncproto.h
UVES_combsynthThAr.o: /opt/local/include/X11/Xosdefs.h
UVES_combsynthThAr.o: /opt/local/include/X11/Xutil.h
UVES_combsynthThAr.o: /opt/local/include/X11/keysym.h
UVES_combsynthThAr.o: /opt/local/include/X11/keysymdef.h
UVES_combsynthThAr.o: /opt/local/include/X11/Xatom.h
UVES_combsynthThAr.o: /opt/local/include/X11/Xos.h
UVES_combsynthThAr.o: /opt/local/include/X11/Xarch.h
UVES_combsynthThAr.o: /opt/local/include/X11/extensions/shape.h
UVES_combsynthThAr.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_combsynthThAr.o: charstr.h utils.h gamm.h const.h error.h
UVES_confit.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_confit.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_confit.o: /opt/local/include/X11/Xfuncproto.h
UVES_confit.o: /opt/local/include/X11/Xosdefs.h
UVES_confit.o: /opt/local/include/X11/Xutil.h /opt/local/include/X11/keysym.h
UVES_confit.o: /opt/local/include/X11/keysymdef.h
UVES_confit.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_confit.o: /opt/local/include/X11/Xarch.h
UVES_confit.o: /opt/local/include/X11/extensions/shape.h
UVES_confit.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_confit.o: utils.h fit.h sort.h memory.h error.h
UVES_cspec_cont.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_cspec_cont.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_cspec_cont.o: /opt/local/include/X11/Xfuncproto.h
UVES_cspec_cont.o: /opt/local/include/X11/Xosdefs.h
UVES_cspec_cont.o: /opt/local/include/X11/Xutil.h
UVES_cspec_cont.o: /opt/local/include/X11/keysym.h
UVES_cspec_cont.o: /opt/local/include/X11/keysymdef.h
UVES_cspec_cont.o: /opt/local/include/X11/Xatom.h
UVES_cspec_cont.o: /opt/local/include/X11/Xos.h
UVES_cspec_cont.o: /opt/local/include/X11/Xarch.h
UVES_cspec_cont.o: /opt/local/include/X11/extensions/shape.h
UVES_cspec_cont.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_cspec_cont.o: utils.h stats.h const.h error.h
UVES_hex2dec.o: error.h
UVES_hirx_blzfit.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_hirx_blzfit.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_hirx_blzfit.o: /opt/local/include/X11/Xfuncproto.h
UVES_hirx_blzfit.o: /opt/local/include/X11/Xosdefs.h
UVES_hirx_blzfit.o: /opt/local/include/X11/Xutil.h
UVES_hirx_blzfit.o: /opt/local/include/X11/keysym.h
UVES_hirx_blzfit.o: /opt/local/include/X11/keysymdef.h
UVES_hirx_blzfit.o: /opt/local/include/X11/Xatom.h
UVES_hirx_blzfit.o: /opt/local/include/X11/Xos.h
UVES_hirx_blzfit.o: /opt/local/include/X11/Xarch.h
UVES_hirx_blzfit.o: /opt/local/include/X11/extensions/shape.h
UVES_hirx_blzfit.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_hirx_blzfit.o: utils.h fit.h stats.h memory.h error.h
UVES_init_cspec.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_init_cspec.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_init_cspec.o: /opt/local/include/X11/Xfuncproto.h
UVES_init_cspec.o: /opt/local/include/X11/Xosdefs.h
UVES_init_cspec.o: /opt/local/include/X11/Xutil.h
UVES_init_cspec.o: /opt/local/include/X11/keysym.h
UVES_init_cspec.o: /opt/local/include/X11/keysymdef.h
UVES_init_cspec.o: /opt/local/include/X11/Xatom.h
UVES_init_cspec.o: /opt/local/include/X11/Xos.h
UVES_init_cspec.o: /opt/local/include/X11/Xarch.h
UVES_init_cspec.o: /opt/local/include/X11/extensions/shape.h
UVES_init_cspec.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_init_cspec.o: utils.h const.h memory.h error.h
UVES_memspec.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_memspec.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_memspec.o: /opt/local/include/X11/Xfuncproto.h
UVES_memspec.o: /opt/local/include/X11/Xosdefs.h
UVES_memspec.o: /opt/local/include/X11/Xutil.h
UVES_memspec.o: /opt/local/include/X11/keysym.h
UVES_memspec.o: /opt/local/include/X11/keysymdef.h
UVES_memspec.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_memspec.o: /opt/local/include/X11/Xarch.h
UVES_memspec.o: /opt/local/include/X11/extensions/shape.h
UVES_memspec.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_memspec.o: utils.h memory.h error.h
UVES_merge_thar.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_merge_thar.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_merge_thar.o: /opt/local/include/X11/Xfuncproto.h
UVES_merge_thar.o: /opt/local/include/X11/Xosdefs.h
UVES_merge_thar.o: /opt/local/include/X11/Xutil.h
UVES_merge_thar.o: /opt/local/include/X11/keysym.h
UVES_merge_thar.o: /opt/local/include/X11/keysymdef.h
UVES_merge_thar.o: /opt/local/include/X11/Xatom.h
UVES_merge_thar.o: /opt/local/include/X11/Xos.h
UVES_merge_thar.o: /opt/local/include/X11/Xarch.h
UVES_merge_thar.o: /opt/local/include/X11/extensions/shape.h
UVES_merge_thar.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_merge_thar.o: utils.h const.h memory.h error.h
UVES_model_resol.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_model_resol.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_model_resol.o: /opt/local/include/X11/Xfuncproto.h
UVES_model_resol.o: /opt/local/include/X11/Xosdefs.h
UVES_model_resol.o: /opt/local/include/X11/Xutil.h
UVES_model_resol.o: /opt/local/include/X11/keysym.h
UVES_model_resol.o: /opt/local/include/X11/keysymdef.h
UVES_model_resol.o: /opt/local/include/X11/Xatom.h
UVES_model_resol.o: /opt/local/include/X11/Xos.h
UVES_model_resol.o: /opt/local/include/X11/Xarch.h
UVES_model_resol.o: /opt/local/include/X11/extensions/shape.h
UVES_model_resol.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_model_resol.o: utils.h fit.h stats.h memory.h const.h error.h
UVES_order_cont.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_order_cont.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_order_cont.o: /opt/local/include/X11/Xfuncproto.h
UVES_order_cont.o: /opt/local/include/X11/Xosdefs.h
UVES_order_cont.o: /opt/local/include/X11/Xutil.h
UVES_order_cont.o: /opt/local/include/X11/keysym.h
UVES_order_cont.o: /opt/local/include/X11/keysymdef.h
UVES_order_cont.o: /opt/local/include/X11/Xatom.h
UVES_order_cont.o: /opt/local/include/X11/Xos.h
UVES_order_cont.o: /opt/local/include/X11/Xarch.h
UVES_order_cont.o: /opt/local/include/X11/extensions/shape.h
UVES_order_cont.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_order_cont.o: utils.h error.h
UVES_order_rejsigedge.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_order_rejsigedge.o: /opt/local/include/X11/Xlib.h
UVES_order_rejsigedge.o: /opt/local/include/X11/X.h
UVES_order_rejsigedge.o: /opt/local/include/X11/Xfuncproto.h
UVES_order_rejsigedge.o: /opt/local/include/X11/Xosdefs.h
UVES_order_rejsigedge.o: /opt/local/include/X11/Xutil.h
UVES_order_rejsigedge.o: /opt/local/include/X11/keysym.h
UVES_order_rejsigedge.o: /opt/local/include/X11/keysymdef.h
UVES_order_rejsigedge.o: /opt/local/include/X11/Xatom.h
UVES_order_rejsigedge.o: /opt/local/include/X11/Xos.h
UVES_order_rejsigedge.o: /opt/local/include/X11/Xarch.h
UVES_order_rejsigedge.o: /opt/local/include/X11/extensions/shape.h
UVES_order_rejsigedge.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_order_rejsigedge.o: charstr.h utils.h memory.h stats.h error.h
UVES_order_sigclip.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_order_sigclip.o: /opt/local/include/X11/Xlib.h
UVES_order_sigclip.o: /opt/local/include/X11/X.h
UVES_order_sigclip.o: /opt/local/include/X11/Xfuncproto.h
UVES_order_sigclip.o: /opt/local/include/X11/Xosdefs.h
UVES_order_sigclip.o: /opt/local/include/X11/Xutil.h
UVES_order_sigclip.o: /opt/local/include/X11/keysym.h
UVES_order_sigclip.o: /opt/local/include/X11/keysymdef.h
UVES_order_sigclip.o: /opt/local/include/X11/Xatom.h
UVES_order_sigclip.o: /opt/local/include/X11/Xos.h
UVES_order_sigclip.o: /opt/local/include/X11/Xarch.h
UVES_order_sigclip.o: /opt/local/include/X11/extensions/shape.h
UVES_order_sigclip.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_order_sigclip.o: charstr.h utils.h stats.h memory.h error.h
UVES_order_stats.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_order_stats.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_order_stats.o: /opt/local/include/X11/Xfuncproto.h
UVES_order_stats.o: /opt/local/include/X11/Xosdefs.h
UVES_order_stats.o: /opt/local/include/X11/Xutil.h
UVES_order_stats.o: /opt/local/include/X11/keysym.h
UVES_order_stats.o: /opt/local/include/X11/keysymdef.h
UVES_order_stats.o: /opt/local/include/X11/Xatom.h
UVES_order_stats.o: /opt/local/include/X11/Xos.h
UVES_order_stats.o: /opt/local/include/X11/Xarch.h
UVES_order_stats.o: /opt/local/include/X11/extensions/shape.h
UVES_order_stats.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_order_stats.o: utils.h stats.h memory.h error.h
UVES_params_init.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_params_init.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_params_init.o: /opt/local/include/X11/Xfuncproto.h
UVES_params_init.o: /opt/local/include/X11/Xosdefs.h
UVES_params_init.o: /opt/local/include/X11/Xutil.h
UVES_params_init.o: /opt/local/include/X11/keysym.h
UVES_params_init.o: /opt/local/include/X11/keysymdef.h
UVES_params_init.o: /opt/local/include/X11/Xatom.h
UVES_params_init.o: /opt/local/include/X11/Xos.h
UVES_params_init.o: /opt/local/include/X11/Xarch.h
UVES_params_init.o: /opt/local/include/X11/extensions/shape.h
UVES_params_init.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_params_init.o: utils.h
UVES_params_set.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_params_set.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_params_set.o: /opt/local/include/X11/Xfuncproto.h
UVES_params_set.o: /opt/local/include/X11/Xosdefs.h
UVES_params_set.o: /opt/local/include/X11/Xutil.h
UVES_params_set.o: /opt/local/include/X11/keysym.h
UVES_params_set.o: /opt/local/include/X11/keysymdef.h
UVES_params_set.o: /opt/local/include/X11/Xatom.h
UVES_params_set.o: /opt/local/include/X11/Xos.h
UVES_params_set.o: /opt/local/include/X11/Xarch.h
UVES_params_set.o: /opt/local/include/X11/extensions/shape.h
UVES_params_set.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_params_set.o: utils.h
UVES_past_actions.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_past_actions.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_past_actions.o: /opt/local/include/X11/Xfuncproto.h
UVES_past_actions.o: /opt/local/include/X11/Xosdefs.h
UVES_past_actions.o: /opt/local/include/X11/Xutil.h
UVES_past_actions.o: /opt/local/include/X11/keysym.h
UVES_past_actions.o: /opt/local/include/X11/keysymdef.h
UVES_past_actions.o: /opt/local/include/X11/Xatom.h
UVES_past_actions.o: /opt/local/include/X11/Xos.h
UVES_past_actions.o: /opt/local/include/X11/Xarch.h
UVES_past_actions.o: /opt/local/include/X11/extensions/shape.h
UVES_past_actions.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_past_actions.o: utils.h fit.h memory.h error.h
UVES_pgenv_init.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_pgenv_init.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_pgenv_init.o: /opt/local/include/X11/Xfuncproto.h
UVES_pgenv_init.o: /opt/local/include/X11/Xosdefs.h
UVES_pgenv_init.o: /opt/local/include/X11/Xutil.h
UVES_pgenv_init.o: /opt/local/include/X11/keysym.h
UVES_pgenv_init.o: /opt/local/include/X11/keysymdef.h
UVES_pgenv_init.o: /opt/local/include/X11/Xatom.h
UVES_pgenv_init.o: /opt/local/include/X11/Xos.h
UVES_pgenv_init.o: /opt/local/include/X11/Xarch.h
UVES_pgenv_init.o: /opt/local/include/X11/extensions/shape.h
UVES_pgenv_init.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_pgenv_init.o: utils.h
UVES_pixscal.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_pixscal.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_pixscal.o: /opt/local/include/X11/Xfuncproto.h
UVES_pixscal.o: /opt/local/include/X11/Xosdefs.h
UVES_pixscal.o: /opt/local/include/X11/Xutil.h
UVES_pixscal.o: /opt/local/include/X11/keysym.h
UVES_pixscal.o: /opt/local/include/X11/keysymdef.h
UVES_pixscal.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_pixscal.o: /opt/local/include/X11/Xarch.h
UVES_pixscal.o: /opt/local/include/X11/extensions/shape.h
UVES_pixscal.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_pixscal.o: utils.h const.h
UVES_plot_cspec.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_plot_cspec.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_plot_cspec.o: /opt/local/include/X11/Xfuncproto.h
UVES_plot_cspec.o: /opt/local/include/X11/Xosdefs.h
UVES_plot_cspec.o: /opt/local/include/X11/Xutil.h
UVES_plot_cspec.o: /opt/local/include/X11/keysym.h
UVES_plot_cspec.o: /opt/local/include/X11/keysymdef.h
UVES_plot_cspec.o: /opt/local/include/X11/Xatom.h
UVES_plot_cspec.o: /opt/local/include/X11/Xos.h
UVES_plot_cspec.o: /opt/local/include/X11/Xarch.h
UVES_plot_cspec.o: /opt/local/include/X11/extensions/shape.h
UVES_plot_cspec.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_plot_cspec.o: utils.h fit.h stats.h sort.h gamm.h input.h memory.h
UVES_plot_cspec.o: const.h error.h
UVES_plot_replay.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_plot_replay.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_plot_replay.o: /opt/local/include/X11/Xfuncproto.h
UVES_plot_replay.o: /opt/local/include/X11/Xosdefs.h
UVES_plot_replay.o: /opt/local/include/X11/Xutil.h
UVES_plot_replay.o: /opt/local/include/X11/keysym.h
UVES_plot_replay.o: /opt/local/include/X11/keysymdef.h
UVES_plot_replay.o: /opt/local/include/X11/Xatom.h
UVES_plot_replay.o: /opt/local/include/X11/Xos.h
UVES_plot_replay.o: /opt/local/include/X11/Xarch.h
UVES_plot_replay.o: /opt/local/include/X11/extensions/shape.h
UVES_plot_replay.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_plot_replay.o: utils.h input.h error.h
UVES_popler.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_popler.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_popler.o: /opt/local/include/X11/Xfuncproto.h
UVES_popler.o: /opt/local/include/X11/Xosdefs.h
UVES_popler.o: /opt/local/include/X11/Xutil.h /opt/local/include/X11/keysym.h
UVES_popler.o: /opt/local/include/X11/keysymdef.h
UVES_popler.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_popler.o: /opt/local/include/X11/Xarch.h
UVES_popler.o: /opt/local/include/X11/extensions/shape.h
UVES_popler.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_popler.o: utils.h stats.h input.h error.h
UVES_r1Dspec.o: /opt/local/include/fitsio.h /opt/local/include/longnam.h
UVES_r1Dspec.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_r1Dspec.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_r1Dspec.o: /opt/local/include/X11/Xfuncproto.h
UVES_r1Dspec.o: /opt/local/include/X11/Xosdefs.h
UVES_r1Dspec.o: /opt/local/include/X11/Xutil.h
UVES_r1Dspec.o: /opt/local/include/X11/keysym.h
UVES_r1Dspec.o: /opt/local/include/X11/keysymdef.h
UVES_r1Dspec.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_r1Dspec.o: /opt/local/include/X11/Xarch.h
UVES_r1Dspec.o: /opt/local/include/X11/extensions/shape.h
UVES_r1Dspec.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_r1Dspec.o: utils.h const.h file.h memory.h error.h
UVES_r2Dspec.o: /opt/local/include/fitsio.h /opt/local/include/longnam.h
UVES_r2Dspec.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_r2Dspec.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_r2Dspec.o: /opt/local/include/X11/Xfuncproto.h
UVES_r2Dspec.o: /opt/local/include/X11/Xosdefs.h
UVES_r2Dspec.o: /opt/local/include/X11/Xutil.h
UVES_r2Dspec.o: /opt/local/include/X11/keysym.h
UVES_r2Dspec.o: /opt/local/include/X11/keysymdef.h
UVES_r2Dspec.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_r2Dspec.o: /opt/local/include/X11/Xarch.h
UVES_r2Dspec.o: /opt/local/include/X11/extensions/shape.h
UVES_r2Dspec.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_r2Dspec.o: utils.h stats.h memory.h error.h const.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/fitsio.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/longnam.h UVES_popler.h pg_plot.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/cpgplot.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/Xlib.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/X.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/Xfuncproto.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/Xosdefs.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/Xutil.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/keysym.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/keysymdef.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/Xatom.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/Xos.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/Xarch.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/extensions/shape.h
UVES_r2Dspec_ESOmer.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_r2Dspec_ESOmer.o: charstr.h utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_harps.o: /opt/local/include/fitsio.h
UVES_r2Dspec_harps.o: /opt/local/include/longnam.h UVES_popler.h pg_plot.h
UVES_r2Dspec_harps.o: /opt/local/include/cpgplot.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/Xlib.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/X.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/Xfuncproto.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/Xosdefs.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/Xutil.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/keysym.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/keysymdef.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/Xatom.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/Xos.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/Xarch.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/extensions/shape.h
UVES_r2Dspec_harps.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_r2Dspec_harps.o: charstr.h utils.h stats.h memory.h error.h const.h
UVES_r2Dspec_hirx.o: /opt/local/include/fitsio.h /opt/local/include/longnam.h
UVES_r2Dspec_hirx.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/Xfuncproto.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/Xosdefs.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/Xutil.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/keysym.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/keysymdef.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/Xatom.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/Xos.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/Xarch.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/extensions/shape.h
UVES_r2Dspec_hirx.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_r2Dspec_hirx.o: utils.h fit.h memory.h error.h const.h
UVES_r2Dspec_iraf.o: /opt/local/include/fitsio.h /opt/local/include/longnam.h
UVES_r2Dspec_iraf.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/Xfuncproto.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/Xosdefs.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/Xutil.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/keysym.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/keysymdef.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/Xatom.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/Xos.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/Xarch.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/extensions/shape.h
UVES_r2Dspec_iraf.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_r2Dspec_iraf.o: utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_iresi.o: /opt/local/include/fitsio.h
UVES_r2Dspec_iresi.o: /opt/local/include/longnam.h UVES_popler.h pg_plot.h
UVES_r2Dspec_iresi.o: /opt/local/include/cpgplot.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/Xlib.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/X.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/Xfuncproto.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/Xosdefs.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/Xutil.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/keysym.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/keysymdef.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/Xatom.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/Xos.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/Xarch.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/extensions/shape.h
UVES_r2Dspec_iresi.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_r2Dspec_iresi.o: charstr.h utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_irls.o: /opt/local/include/fitsio.h /opt/local/include/longnam.h
UVES_r2Dspec_irls.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/Xfuncproto.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/Xosdefs.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/Xutil.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/keysym.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/keysymdef.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/Xatom.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/Xos.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/Xarch.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/extensions/shape.h
UVES_r2Dspec_irls.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_r2Dspec_irls.o: utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_mage.o: /opt/local/include/fitsio.h /opt/local/include/longnam.h
UVES_r2Dspec_mage.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/Xfuncproto.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/Xosdefs.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/Xutil.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/keysym.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/keysymdef.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/Xatom.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/Xos.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/Xarch.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/extensions/shape.h
UVES_r2Dspec_mage.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_r2Dspec_mage.o: utils.h astron.h const.h memory.h error.h
UVES_r2Dspec_makee.o: /opt/local/include/fitsio.h
UVES_r2Dspec_makee.o: /opt/local/include/longnam.h UVES_popler.h pg_plot.h
UVES_r2Dspec_makee.o: /opt/local/include/cpgplot.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/Xlib.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/X.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/Xfuncproto.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/Xosdefs.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/Xutil.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/keysym.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/keysymdef.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/Xatom.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/Xos.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/Xarch.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/extensions/shape.h
UVES_r2Dspec_makee.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_r2Dspec_makee.o: charstr.h utils.h memory.h const.h error.h
UVES_ratmask.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_ratmask.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_ratmask.o: /opt/local/include/X11/Xfuncproto.h
UVES_ratmask.o: /opt/local/include/X11/Xosdefs.h
UVES_ratmask.o: /opt/local/include/X11/Xutil.h
UVES_ratmask.o: /opt/local/include/X11/keysym.h
UVES_ratmask.o: /opt/local/include/X11/keysymdef.h
UVES_ratmask.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_ratmask.o: /opt/local/include/X11/Xarch.h
UVES_ratmask.o: /opt/local/include/X11/extensions/shape.h
UVES_ratmask.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_ratmask.o: utils.h memory.h file.h error.h
UVES_redispers.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_redispers.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_redispers.o: /opt/local/include/X11/Xfuncproto.h
UVES_redispers.o: /opt/local/include/X11/Xosdefs.h
UVES_redispers.o: /opt/local/include/X11/Xutil.h
UVES_redispers.o: /opt/local/include/X11/keysym.h
UVES_redispers.o: /opt/local/include/X11/keysymdef.h
UVES_redispers.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_redispers.o: /opt/local/include/X11/Xarch.h
UVES_redispers.o: /opt/local/include/X11/extensions/shape.h
UVES_redispers.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_redispers.o: utils.h memory.h error.h
UVES_replay_control.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_replay_control.o: /opt/local/include/X11/Xlib.h
UVES_replay_control.o: /opt/local/include/X11/X.h
UVES_replay_control.o: /opt/local/include/X11/Xfuncproto.h
UVES_replay_control.o: /opt/local/include/X11/Xosdefs.h
UVES_replay_control.o: /opt/local/include/X11/Xutil.h
UVES_replay_control.o: /opt/local/include/X11/keysym.h
UVES_replay_control.o: /opt/local/include/X11/keysymdef.h
UVES_replay_control.o: /opt/local/include/X11/Xatom.h
UVES_replay_control.o: /opt/local/include/X11/Xos.h
UVES_replay_control.o: /opt/local/include/X11/Xarch.h
UVES_replay_control.o: /opt/local/include/X11/extensions/shape.h
UVES_replay_control.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_replay_control.o: charstr.h utils.h input.h stats.h memory.h error.h
UVES_replay_control.o: const.h
UVES_rescale_region.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_rescale_region.o: /opt/local/include/X11/Xlib.h
UVES_rescale_region.o: /opt/local/include/X11/X.h
UVES_rescale_region.o: /opt/local/include/X11/Xfuncproto.h
UVES_rescale_region.o: /opt/local/include/X11/Xosdefs.h
UVES_rescale_region.o: /opt/local/include/X11/Xutil.h
UVES_rescale_region.o: /opt/local/include/X11/keysym.h
UVES_rescale_region.o: /opt/local/include/X11/keysymdef.h
UVES_rescale_region.o: /opt/local/include/X11/Xatom.h
UVES_rescale_region.o: /opt/local/include/X11/Xos.h
UVES_rescale_region.o: /opt/local/include/X11/Xarch.h
UVES_rescale_region.o: /opt/local/include/X11/extensions/shape.h
UVES_rescale_region.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_rescale_region.o: charstr.h utils.h stats.h fit.h memory.h error.h
UVES_revwpol.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_revwpol.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_revwpol.o: /opt/local/include/X11/Xfuncproto.h
UVES_revwpol.o: /opt/local/include/X11/Xosdefs.h
UVES_revwpol.o: /opt/local/include/X11/Xutil.h
UVES_revwpol.o: /opt/local/include/X11/keysym.h
UVES_revwpol.o: /opt/local/include/X11/keysymdef.h
UVES_revwpol.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_revwpol.o: /opt/local/include/X11/Xarch.h
UVES_revwpol.o: /opt/local/include/X11/extensions/shape.h
UVES_revwpol.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_revwpol.o: utils.h fit.h memory.h error.h
UVES_rinputfile.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_rinputfile.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_rinputfile.o: /opt/local/include/X11/Xfuncproto.h
UVES_rinputfile.o: /opt/local/include/X11/Xosdefs.h
UVES_rinputfile.o: /opt/local/include/X11/Xutil.h
UVES_rinputfile.o: /opt/local/include/X11/keysym.h
UVES_rinputfile.o: /opt/local/include/X11/keysymdef.h
UVES_rinputfile.o: /opt/local/include/X11/Xatom.h
UVES_rinputfile.o: /opt/local/include/X11/Xos.h
UVES_rinputfile.o: /opt/local/include/X11/Xarch.h
UVES_rinputfile.o: /opt/local/include/X11/extensions/shape.h
UVES_rinputfile.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_rinputfile.o: utils.h file.h error.h
UVES_rFITSlist.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_rFITSlist.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_rFITSlist.o: /opt/local/include/X11/Xfuncproto.h
UVES_rFITSlist.o: /opt/local/include/X11/Xosdefs.h
UVES_rFITSlist.o: /opt/local/include/X11/Xutil.h
UVES_rFITSlist.o: /opt/local/include/X11/keysym.h
UVES_rFITSlist.o: /opt/local/include/X11/keysymdef.h
UVES_rFITSlist.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_rFITSlist.o: /opt/local/include/X11/Xarch.h
UVES_rFITSlist.o: /opt/local/include/X11/extensions/shape.h
UVES_rFITSlist.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_rFITSlist.o: utils.h file.h error.h
UVES_rMacmap.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_rMacmap.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_rMacmap.o: /opt/local/include/X11/Xfuncproto.h
UVES_rMacmap.o: /opt/local/include/X11/Xosdefs.h
UVES_rMacmap.o: /opt/local/include/X11/Xutil.h
UVES_rMacmap.o: /opt/local/include/X11/keysym.h
UVES_rMacmap.o: /opt/local/include/X11/keysymdef.h
UVES_rMacmap.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_rMacmap.o: /opt/local/include/X11/Xarch.h
UVES_rMacmap.o: /opt/local/include/X11/extensions/shape.h
UVES_rMacmap.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_rMacmap.o: utils.h memory.h file.h error.h
UVES_rUPLfile.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_rUPLfile.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_rUPLfile.o: /opt/local/include/X11/Xfuncproto.h
UVES_rUPLfile.o: /opt/local/include/X11/Xosdefs.h
UVES_rUPLfile.o: /opt/local/include/X11/Xutil.h
UVES_rUPLfile.o: /opt/local/include/X11/keysym.h
UVES_rUPLfile.o: /opt/local/include/X11/keysymdef.h
UVES_rUPLfile.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_rUPLfile.o: /opt/local/include/X11/Xarch.h
UVES_rUPLfile.o: /opt/local/include/X11/extensions/shape.h
UVES_rUPLfile.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_rUPLfile.o: utils.h file.h error.h
UVES_rvshift.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_rvshift.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_rvshift.o: /opt/local/include/X11/Xfuncproto.h
UVES_rvshift.o: /opt/local/include/X11/Xosdefs.h
UVES_rvshift.o: /opt/local/include/X11/Xutil.h
UVES_rvshift.o: /opt/local/include/X11/keysym.h
UVES_rvshift.o: /opt/local/include/X11/keysymdef.h
UVES_rvshift.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_rvshift.o: /opt/local/include/X11/Xarch.h
UVES_rvshift.o: /opt/local/include/X11/extensions/shape.h
UVES_rvshift.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_rvshift.o: utils.h file.h error.h
UVES_select_subspec.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_select_subspec.o: /opt/local/include/X11/Xlib.h
UVES_select_subspec.o: /opt/local/include/X11/X.h
UVES_select_subspec.o: /opt/local/include/X11/Xfuncproto.h
UVES_select_subspec.o: /opt/local/include/X11/Xosdefs.h
UVES_select_subspec.o: /opt/local/include/X11/Xutil.h
UVES_select_subspec.o: /opt/local/include/X11/keysym.h
UVES_select_subspec.o: /opt/local/include/X11/keysymdef.h
UVES_select_subspec.o: /opt/local/include/X11/Xatom.h
UVES_select_subspec.o: /opt/local/include/X11/Xos.h
UVES_select_subspec.o: /opt/local/include/X11/Xarch.h
UVES_select_subspec.o: /opt/local/include/X11/extensions/shape.h
UVES_select_subspec.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_select_subspec.o: charstr.h utils.h memory.h error.h
UVES_set_wavelen_scale.o: UVES_popler.h pg_plot.h
UVES_set_wavelen_scale.o: /opt/local/include/cpgplot.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/Xlib.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/X.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/Xfuncproto.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/Xosdefs.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/Xutil.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/keysym.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/keysymdef.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/Xatom.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/Xos.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/Xarch.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/extensions/shape.h
UVES_set_wavelen_scale.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_set_wavelen_scale.o: charstr.h utils.h const.h memory.h error.h
UVES_skysub.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_skysub.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_skysub.o: /opt/local/include/X11/Xfuncproto.h
UVES_skysub.o: /opt/local/include/X11/Xosdefs.h
UVES_skysub.o: /opt/local/include/X11/Xutil.h /opt/local/include/X11/keysym.h
UVES_skysub.o: /opt/local/include/X11/keysymdef.h
UVES_skysub.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_skysub.o: /opt/local/include/X11/Xarch.h
UVES_skysub.o: /opt/local/include/X11/extensions/shape.h
UVES_skysub.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_skysub.o: utils.h memory.h error.h
UVES_synthThAr.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_synthThAr.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_synthThAr.o: /opt/local/include/X11/Xfuncproto.h
UVES_synthThAr.o: /opt/local/include/X11/Xosdefs.h
UVES_synthThAr.o: /opt/local/include/X11/Xutil.h
UVES_synthThAr.o: /opt/local/include/X11/keysym.h
UVES_synthThAr.o: /opt/local/include/X11/keysymdef.h
UVES_synthThAr.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_synthThAr.o: /opt/local/include/X11/Xarch.h
UVES_synthThAr.o: /opt/local/include/X11/extensions/shape.h
UVES_synthThAr.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_synthThAr.o: utils.h gamm.h const.h error.h
UVES_thar_sigarray.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_thar_sigarray.o: /opt/local/include/X11/Xlib.h
UVES_thar_sigarray.o: /opt/local/include/X11/X.h
UVES_thar_sigarray.o: /opt/local/include/X11/Xfuncproto.h
UVES_thar_sigarray.o: /opt/local/include/X11/Xosdefs.h
UVES_thar_sigarray.o: /opt/local/include/X11/Xutil.h
UVES_thar_sigarray.o: /opt/local/include/X11/keysym.h
UVES_thar_sigarray.o: /opt/local/include/X11/keysymdef.h
UVES_thar_sigarray.o: /opt/local/include/X11/Xatom.h
UVES_thar_sigarray.o: /opt/local/include/X11/Xos.h
UVES_thar_sigarray.o: /opt/local/include/X11/Xarch.h
UVES_thar_sigarray.o: /opt/local/include/X11/extensions/shape.h
UVES_thar_sigarray.o: /opt/local/include/X11/extensions/shapeconst.h
UVES_thar_sigarray.o: charstr.h utils.h memory.h error.h
UVES_undo_lastact.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_undo_lastact.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_undo_lastact.o: /opt/local/include/X11/Xfuncproto.h
UVES_undo_lastact.o: /opt/local/include/X11/Xosdefs.h
UVES_undo_lastact.o: /opt/local/include/X11/Xutil.h
UVES_undo_lastact.o: /opt/local/include/X11/keysym.h
UVES_undo_lastact.o: /opt/local/include/X11/keysymdef.h
UVES_undo_lastact.o: /opt/local/include/X11/Xatom.h
UVES_undo_lastact.o: /opt/local/include/X11/Xos.h
UVES_undo_lastact.o: /opt/local/include/X11/Xarch.h
UVES_undo_lastact.o: /opt/local/include/X11/extensions/shape.h
UVES_undo_lastact.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_undo_lastact.o: utils.h error.h
UVES_vhelio.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_vhelio.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_vhelio.o: /opt/local/include/X11/Xfuncproto.h
UVES_vhelio.o: /opt/local/include/X11/Xosdefs.h
UVES_vhelio.o: /opt/local/include/X11/Xutil.h /opt/local/include/X11/keysym.h
UVES_vhelio.o: /opt/local/include/X11/keysymdef.h
UVES_vhelio.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_vhelio.o: /opt/local/include/X11/Xarch.h
UVES_vhelio.o: /opt/local/include/X11/extensions/shape.h
UVES_vhelio.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_vhelio.o: utils.h astron.h const.h error.h
UVES_wDATfile.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_wDATfile.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_wDATfile.o: /opt/local/include/X11/Xfuncproto.h
UVES_wDATfile.o: /opt/local/include/X11/Xosdefs.h
UVES_wDATfile.o: /opt/local/include/X11/Xutil.h
UVES_wDATfile.o: /opt/local/include/X11/keysym.h
UVES_wDATfile.o: /opt/local/include/X11/keysymdef.h
UVES_wDATfile.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_wDATfile.o: /opt/local/include/X11/Xarch.h
UVES_wDATfile.o: /opt/local/include/X11/extensions/shape.h
UVES_wDATfile.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_wDATfile.o: utils.h file.h error.h
UVES_wFITSfile.o: /opt/local/include/fitsio.h /opt/local/include/longnam.h
UVES_wFITSfile.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_wFITSfile.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_wFITSfile.o: /opt/local/include/X11/Xfuncproto.h
UVES_wFITSfile.o: /opt/local/include/X11/Xosdefs.h
UVES_wFITSfile.o: /opt/local/include/X11/Xutil.h
UVES_wFITSfile.o: /opt/local/include/X11/keysym.h
UVES_wFITSfile.o: /opt/local/include/X11/keysymdef.h
UVES_wFITSfile.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_wFITSfile.o: /opt/local/include/X11/Xarch.h
UVES_wFITSfile.o: /opt/local/include/X11/extensions/shape.h
UVES_wFITSfile.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_wFITSfile.o: utils.h memory.h error.h
UVES_wrawFITS.o: /opt/local/include/fitsio.h /opt/local/include/longnam.h
UVES_wrawFITS.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_wrawFITS.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_wrawFITS.o: /opt/local/include/X11/Xfuncproto.h
UVES_wrawFITS.o: /opt/local/include/X11/Xosdefs.h
UVES_wrawFITS.o: /opt/local/include/X11/Xutil.h
UVES_wrawFITS.o: /opt/local/include/X11/keysym.h
UVES_wrawFITS.o: /opt/local/include/X11/keysymdef.h
UVES_wrawFITS.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_wrawFITS.o: /opt/local/include/X11/Xarch.h
UVES_wrawFITS.o: /opt/local/include/X11/extensions/shape.h
UVES_wrawFITS.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_wrawFITS.o: utils.h memory.h error.h
UVES_wUPLfile.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_wUPLfile.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_wUPLfile.o: /opt/local/include/X11/Xfuncproto.h
UVES_wUPLfile.o: /opt/local/include/X11/Xosdefs.h
UVES_wUPLfile.o: /opt/local/include/X11/Xutil.h
UVES_wUPLfile.o: /opt/local/include/X11/keysym.h
UVES_wUPLfile.o: /opt/local/include/X11/keysymdef.h
UVES_wUPLfile.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_wUPLfile.o: /opt/local/include/X11/Xarch.h
UVES_wUPLfile.o: /opt/local/include/X11/extensions/shape.h
UVES_wUPLfile.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_wUPLfile.o: utils.h file.h error.h
UVES_wpol.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_wpol.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_wpol.o: /opt/local/include/X11/Xfuncproto.h
UVES_wpol.o: /opt/local/include/X11/Xosdefs.h /opt/local/include/X11/Xutil.h
UVES_wpol.o: /opt/local/include/X11/keysym.h
UVES_wpol.o: /opt/local/include/X11/keysymdef.h
UVES_wpol.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_wpol.o: /opt/local/include/X11/Xarch.h
UVES_wpol.o: /opt/local/include/X11/extensions/shape.h
UVES_wpol.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h utils.h
UVES_wpol.o: fit.h stats.h memory.h const.h error.h
UVES_tmpfit.o: UVES_popler.h pg_plot.h /opt/local/include/cpgplot.h
UVES_tmpfit.o: /opt/local/include/X11/Xlib.h /opt/local/include/X11/X.h
UVES_tmpfit.o: /opt/local/include/X11/Xfuncproto.h
UVES_tmpfit.o: /opt/local/include/X11/Xosdefs.h
UVES_tmpfit.o: /opt/local/include/X11/Xutil.h /opt/local/include/X11/keysym.h
UVES_tmpfit.o: /opt/local/include/X11/keysymdef.h
UVES_tmpfit.o: /opt/local/include/X11/Xatom.h /opt/local/include/X11/Xos.h
UVES_tmpfit.o: /opt/local/include/X11/Xarch.h
UVES_tmpfit.o: /opt/local/include/X11/extensions/shape.h
UVES_tmpfit.o: /opt/local/include/X11/extensions/shapeconst.h charstr.h
UVES_tmpfit.o: utils.h fit.h memory.h const.h error.h
zbrent.o: fit.h error.h
