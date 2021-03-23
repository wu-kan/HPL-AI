AC_PREREQ([2.69])

AC_INIT(hplai, 2.3b, i@wu-kan.cn)
AC_CONFIG_SRCDIR([include/hplai.hh])
AC_CONFIG_HEADERS([include/hplai_config.hh])

AX_PROG_CC_MPI

AC_PROG_RANLIB

AC_PROG_INSTALL

AM_INIT_AUTOMAKE([subdir-objects])

AM_PROG_CC_C_O

dnl
dnl AX_BLAS requires Fortran compiler and detects fortran libraries in $FLIBS
dnl
dnl AX_BLAS(LIBS="$BLAS_LIBS $LIBS $FLIBS")
dnl

AC_DEFUN([HPL_BLAS], [

AC_PREREQ(2.69)

hpl_blas_ok=no

dnl FIXME: add --with-blas="<library spec>"

current_LIBS="$LIBS"

cat <<HPLEOF > hplvars.txt
name1=OpenBLAS
rout1=dgemm_
libs1=-lopenblas -lm

name2=Atlas Fortran BLAS
rout2=dgemm_
libs2=-lf77blas -latlas

name3=Sequential Intel MKL LP64 (group)
rout3=dgemm_
libs3=-Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread

name4=Sequential Intel MKL LP64
rout4=dgemm_
libs4=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

name5=AMD's ACML
rout5=dgemm_
libs5=-lacml -lm

name6=Accelerate
rout6=dgemm_
libs6=-framework Accelerate

name7=Apple VecLib
rout7=dgemm_
libs7=-framework vecLib

name8=IBM ESSL
rout8=dgemm_
libs8=-lessl

name9=NVIDIA nvblas
rout9=dgemm_
libs9=-lnvblas

name10=Generic BLAS
rout10=dgemm_
libs10=-lblas

HPLEOF
for hpl_i in 1 2 3 4 5 6 7 8 9 10;
do
if test  x$hpl_blas_ok = xno; then
  name="`grep ^name${hpl_i}= hplvars.txt | sed s/^name${hpl_i}=//`"
  rout="`grep ^rout${hpl_i}= hplvars.txt | sed s/^rout${hpl_i}=//`"
  libs="`grep ^libs${hpl_i}= hplvars.txt | sed s/^libs${hpl_i}=//`"
  AC_MSG_CHECKING([for [$]rout in [$]name])

  LIBS="[$]libs"
  AC_TRY_LINK_FUNC([$]rout, [hpl_blas_ok=yes;BLAS_LIBS="[$]libs"])
  LIBS="$current_LIBS"

  AC_MSG_RESULT($hpl_blas_ok)
fi
done
rm hplvars.txt

if test  x$hpl_blas_ok = xno; then
dnl
AC_MSG_CHECKING([for dgemm_ in OpenBLAS])
AC_CHECK_LIB(openblas, dgemm_, [hpl_blas_ok=yes;BLAS_LIBS="-lopenblas"])
AC_MSG_RESULT($hpl_blas_ok)
dnl
fi

AC_SUBST(BLAS_LIBS)

# If present, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$hpl_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        hpl_blas_ok=no
        $2
fi

])dnl HPL_BLAS

HPL_BLAS(LIBS="$BLAS_LIBS $LIBS",AC_MSG_ERROR([BLAS not found]))

dnl FIXME: test for CBLAS: Atlas, MKL, OpenBLAS, ESSL, ...
dnl FIXME: test for GSL CBLAS

AC_CHECK_FUNCS([dgemm_])

AC_CHECK_HEADERS([mpi.h])

AC_CHECK_HEADERS([blas.hh])

AC_PROG_CXX

dnl FIXME: AX_CXX_CHECK_LIB(blaspp, [blas::blaspp_version()])
AC_CHECK_LIB(blaspp, dgemm_)

AC_CONFIG_FILES([Makefile src/Makefile testing/Makefile])

AC_OUTPUT
