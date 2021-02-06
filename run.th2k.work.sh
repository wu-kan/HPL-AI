make clean
aclocal
autoconf
autoheader
automake --add-missing
./configure \
    LIBS="-lmpi" \
    CPPFLAGS="-DHPL_CALL_CBLAS -DHPL_DETAILED_TIMING"
make -j
cd testing/ptest
mpiexec -n 4 ../xhplai