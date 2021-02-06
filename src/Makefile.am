AM_CPPFLAGS = -I$(top_srcdir)/../include

lib_LIBRARIES = libhpl.a libhplai.a

libhplai_a_SOURCES = \
auxil/HPLAI_alatcpy.cpp auxil/HPLAI_alacpy.cpp auxil/HPLAI_alamch.cpp \
blas/HPLAI_blas.cpp \
blas/HPLAI_ascal.cpp blas/HPLAI_atrsm.cpp blas/HPLAI_atrsv.cpp blas/HPLAI_iaamax.cpp \
blas/HPLAI_agemv.cpp blas/HPLAI_aaxpy.cpp \
blas/HPLAI_acopy.cpp blas/HPLAI_agemm.cpp blas/HPLAI_ager.cpp \
comm/HPLAI_sdrv.cpp comm/HPLAI_send.cpp comm/HPLAI_recv.cpp comm/HPLAI_bcast.cpp \
comm/HPLAI_binit.cpp comm/HPLAI_bwait.cpp comm/HPLAI_blong.cpp comm/HPLAI_1ring.cpp \
comm/HPLAI_1rinM.cpp comm/HPLAI_2rinM.cpp comm/HPLAI_2ring.cpp comm/HPLAI_blonM.cpp comm/HPLAI_packL.cpp \
comm/HPLAI_copyL.cpp \
grid/HPLAI_broadcast_afloat.cpp \
panel/HPLAI_papanel_disp.cpp panel/HPLAI_papanel_free.cpp panel/HPLAI_papanel_init.cpp panel/HPLAI_papanel_new.cpp \
pauxil/HPL_indxl2g.c pauxil/HPL_indxg2lp.c \
pfact/HPLAI_pafact.cpp \
pgesv/HPLAI_pdgesv.cpp pgesv/HPLAI_pagesv.cpp pgesv/HPLAI_patrsv.cpp

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
