AM_CPPFLAGS = -I$(top_srcdir)/../include

lib_LIBRARIES = libhpl.a libhplai.a

libhplai_a_SOURCES = \
auxil/HPLAI_alatcpy.cc auxil/HPLAI_alacpy.cc auxil/HPLAI_alamch.cc \
blas/HPLAI_blas.cc \
comm/HPLAI_sdrv.cc comm/HPLAI_send.cc comm/HPLAI_recv.cc comm/HPLAI_bcast.cc \
comm/HPLAI_binit.cc comm/HPLAI_bwait.cc comm/HPLAI_blong.cc comm/HPLAI_1ring.cc \
comm/HPLAI_1rinM.cc comm/HPLAI_2rinM.cc comm/HPLAI_2ring.cc comm/HPLAI_blonM.cc comm/HPLAI_packL.cc \
comm/HPLAI_copyL.cc \
grid/HPLAI_reduce_AFLOAT.cc grid/HPLAI_sum_AFLOAT.cc \
grid/HPLAI_all_reduce_AFLOAT.cc grid/HPLAI_broadcast_AFLOAT.cc grid/HPLAI_max_AFLOAT.cc \
grid/HPLAI_min_AFLOAT.cc \
panel/HPLAI_papanel_disp.cc panel/HPLAI_papanel_free.cc panel/HPLAI_papanel_init.cc panel/HPLAI_papanel_new.cc \
pauxil/HPL_indxl2g.c pauxil/HPL_indxg2lp.c \
pauxil/HPLAI_palamch.cc pauxil/HPLAI_palange.cc \
pauxil/HPLAI_alaswp00N.cc pauxil/HPLAI_alaswp01N.cc pauxil/HPLAI_alaswp01T.cc \
pauxil/HPLAI_alaswp02N.cc pauxil/HPLAI_alaswp03N.cc pauxil/HPLAI_alaswp03T.cc \
pauxil/HPLAI_alaswp04N.cc pauxil/HPLAI_alaswp04T.cc pauxil/HPLAI_alaswp05N.cc \
pauxil/HPLAI_alaswp05T.cc pauxil/HPLAI_alaswp06N.cc pauxil/HPLAI_alaswp06T.cc \
pauxil/HPLAI_alaswp10N.cc \
pfact/HPLAI_papanllN.cc pfact/HPLAI_papanllT.cc pfact/HPLAI_papanrlN.cc \
pfact/HPLAI_papanrlT.cc pfact/HPLAI_parpancrN.cc pfact/HPLAI_parpancrT.cc \
pfact/HPLAI_parpanllN.cc pfact/HPLAI_parpanllT.cc pfact/HPLAI_parpanrlN.cc pfact/HPLAI_parpanrlT.cc \
pfact/HPLAI_pamxswp.cc pfact/HPLAI_pafact.cc pfact/HPLAI_alocmax.cc \
pfact/HPLAI_papancrT.cc pfact/HPLAI_papancrN.cc \
pfact/HPLAI_alocswpN.cc pfact/HPLAI_alocswpT.cc \
pgesv/HPLAI_pdgesv.cc \
pgesv/HPLAI_pagesv0.cc pgesv/HPLAI_pagesv.cc pgesv/HPLAI_pagesvK2.cc \
pgesv/HPLAI_patrsv.cc \
pgesv/HPLAI_paupdateNN.cc pgesv/HPLAI_paupdateNT.cc pgesv/HPLAI_paupdateTN.cc pgesv/HPLAI_paupdateTT.cc \
pgesv/HPLAI_equil.cc pgesv/HPLAI_pipid.cc pgesv/HPLAI_plindx0.cc \
pgesv/HPLAI_plindx10.cc pgesv/HPLAI_plindx1.cc \
pgesv/HPLAI_rollN.cc pgesv/HPLAI_rollT.cc pgesv/HPLAI_spreadN.cc pgesv/HPLAI_spreadT.cc \
pgesv/HPLAI_palaswp00N.cc pgesv/HPLAI_palaswp00T.cc pgesv/HPLAI_palaswp01N.cc pgesv/HPLAI_palaswp01T.cc

libhpl_a_SOURCES = \
auxil/HPL_dlatcpy.c auxil/HPL_fprintf.c auxil/HPL_dlacpy.c auxil/HPL_dlamch.c \
blas/HPL_dscal.c blas/HPL_dtrsm.c blas/HPL_dtrsv.c blas/HPL_idamax.c \
blas/HPL_dgemv.c blas/HPL_dscal.c blas/HPL_daxpy.c \
blas/HPL_dcopy.c blas/HPL_dgemm.c blas/HPL_dgemv.c blas/HPL_dger.c \
comm/HPL_sdrv.c comm/HPL_send.c comm/HPL_recv.c comm/HPL_bcast.c \
comm/HPL_binit.c comm/HPL_bwait.c comm/HPL_blong.c comm/HPL_1ring.c \
comm/HPL_1rinM.c comm/HPL_2rinM.c comm/HPL_2ring.c comm/HPL_blonM.c comm/HPL_packL.c \
grid/HPL_reduce.c grid/HPL_sum.c grid/HPL_grid_info.c grid/HPL_grid_init.c \
grid/HPL_all_reduce.c grid/HPL_broadcast.c grid/HPL_grid_exit.c grid/HPL_max.c \
grid/HPL_min.c grid/HPL_all_reduce.c grid/HPL_barrier.c \
panel/HPL_pdpanel_disp.c panel/HPL_pdpanel_free.c panel/HPL_pdpanel_init.c panel/HPL_pdpanel_new.c \
pauxil/HPL_pdlamch.c pauxil/HPL_pdlange.c \
pauxil/HPL_indxg2p.c pauxil/HPL_numroc.c pauxil/HPL_numrocI.c pauxil/HPL_numrocI.c \
pauxil/HPL_dlaswp00N.c pauxil/HPL_dlaswp01N.c pauxil/HPL_dlaswp01T.c \
pauxil/HPL_dlaswp02N.c pauxil/HPL_dlaswp03N.c pauxil/HPL_dlaswp03T.c \
pauxil/HPL_dlaswp04N.c pauxil/HPL_dlaswp04T.c pauxil/HPL_dlaswp05N.c \
pauxil/HPL_dlaswp05T.c pauxil/HPL_dlaswp06N.c pauxil/HPL_dlaswp06T.c \
pauxil/HPL_infog2l.c pauxil/HPL_dlaswp10N.c pauxil/HPL_pwarn.c \
pfact/HPL_pdpanllN.c pfact/HPL_pdpanllT.c pfact/HPL_pdpanrlN.c \
pfact/HPL_pdpanrlT.c pfact/HPL_pdrpancrN.c pfact/HPL_pdrpancrT.c \
pfact/HPL_pdrpanllN.c pfact/HPL_pdrpanllT.c pfact/HPL_pdrpanrlN.c pfact/HPL_pdrpanrlT.c \
pfact/HPL_pdmxswp.c pfact/HPL_pdfact.c pfact/HPL_dlocmax.c \
pfact/HPL_pdpancrT.c pfact/HPL_pdpancrN.c pfact/HPL_dlocmax.c \
pfact/HPL_dlocswpN.c pfact/HPL_dlocswpT.c pfact/HPL_pdmxswp.c \
pfact/HPL_pdpanllN.c pfact/HPL_pdpanllT.c pfact/HPL_pdpanrlN.c \
pfact/HPL_pdpanrlT.c pfact/HPL_pdrpancrN.c pfact/HPL_pdrpancrT.c \
pfact/HPL_pdrpanllN.c pfact/HPL_pdrpanllT.c pfact/HPL_pdrpanrlN.c \
pfact/HPL_pdrpanrlT.c pauxil/HPL_pabort.c pauxil/HPL_pdlamch.c \
pgesv/HPL_pdgesv0.c pgesv/HPL_pdgesv.c pgesv/HPL_pdgesvK1.c pgesv/HPL_pdgesvK2.c \
pgesv/HPL_pdupdateNN.c pgesv/HPL_pdupdateNT.c pgesv/HPL_pdupdateTN.c pgesv/HPL_pdupdateTT.c \
pgesv/HPL_equil.c pgesv/HPL_pipid.c pgesv/HPL_plindx0.c \
pgesv/HPL_plindx10.c pgesv/HPL_plindx1.c pgesv/HPL_plindx10.c \
pgesv/HPL_rollN.c pgesv/HPL_rollT.c pgesv/HPL_spreadN.c pgesv/HPL_spreadT.c \
pgesv/HPL_pdgesv0.c pgesv/HPL_pdgesv.c pgesv/HPL_pdgesvK1.c pgesv/HPL_pdgesvK2.c pgesv/HPL_pdgesv0.c pgesv/HPL_pdgesvK2.c \
pgesv/HPL_pdlaswp00N.c pgesv/HPL_pdlaswp00T.c pgesv/HPL_pdlaswp01N.c pgesv/HPL_pdlaswp01T.c \
pgesv/HPL_pdtrsv.c pgesv/HPL_pdupdateNN.c pgesv/HPL_pdupdateNT.c pgesv/HPL_pdupdateTN.c \
pgesv/HPL_pdupdateTT.c pgesv/HPL_logsort.c pgesv/HPL_perm.c
