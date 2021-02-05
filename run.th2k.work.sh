make clean
aclocal
autoconf
autoheader
automake
./configure \
    LIBS="-lmpi" \
    CPPFLAGS="-DHPL_CALL_CBLAS -DHPL_DETAILED_TIMING"
make -j
cd testing/ptest
mpiexec -n 4 ../xhpl_ai