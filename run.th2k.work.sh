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
    CPPFLAGS="-DHPL_CALL_CBLAS -DHPL_DETAILED_TIMING -DHPLAI_T_AFLOAT=double"
make -j
echo 'HPLinpack benchmark input file
Innovative Computing Laboratory, University of Tennessee
HPL.out      output file name (if any)
6            device out (6=stdout,7=stderr,file)
1            # of problems sizes (N)
9 29 30 34 35  Ns
1            # of NBs
1 2 3 4      NBs
0            PMAP process mapping (0=Row-,1=Column-major)
1            # of process grids (P x Q)
2 1 4        Ps
2 4 1        Qs
16.0         threshold
1            # of panel fact
0 1 2        PFACTs (0=left, 1=Crout, 2=Right)
1            # of recursive stopping criterium
2 4          NBMINs (>= 1)
1            # of panels in recursion
2            NDIVs
1            # of recursive panel fact.
0 1 2        RFACTs (0=left, 1=Crout, 2=Right)
1            # of broadcast
0            BCASTs (0=1rg,1=1rM,2=2rg,3=2rM,4=Lng,5=LnM)
1            # of lookahead depth
0            DEPTHs (>=0)
2            SWAP (0=bin-exch,1=long,2=mix)
64           swapping threshold
0            L1 in (0=transposed,1=no-transposed) form
0            U  in (0=transposed,1=no-transposed) form
1            Equilibration (0=no,1=yes)
8            memory alignment in double (> 0)' >HPL.dat
mpiexec -n 4 testing/xhpl
mpiexec -n 4 testing/xhplai
git clean -d -f -q -x
