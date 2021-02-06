#!/bin/bash

source ../spack/share/spack/setup-env.sh
spack unload -a
spack load gcc@7.5.0%gcc@4.8.5
spack load autoconf@2.69%gcc@7.5.0
spack load autoconf-archive@2019.01.06%gcc@7.5.0
spack load automake@1.16.3%gcc@7.5.0
spack load openmpi@4.0.5%gcc@7.5.0
spack load openblas@0.3.13%gcc@7.5.0
aclocal
autoconf
autoheader
automake --add-missing
./configure \
    CPPFLAGS="-DHPL_CALL_CBLAS -DHPL_DETAILED_TIMING"
make -j
cd testing/ptest
mpiexec -n 4 ../xhplai
cd ../..
git clean -d -f -q -x