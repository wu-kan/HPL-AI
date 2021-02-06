# HPL-AI

## 开发环境

```bash
source ../spack/share/spack/setup-env.sh
spack unload -a
spack load gcc@7.5.0%gcc@4.8.5
spack load autoconf@2.69%gcc@7.5.0
spack load automake@1.16.3%gcc@7.5.0
spack load openmpi@4.0.5%gcc@7.5.0
spack load openblas@0.3.13%gcc@7.5.0
```

## 发布过程

```bash
aclocal
autoconf
autoheader
automake --add-missing
```

## 编译过程

```bash
./configure \
    LIBS="-lmpi" \
    CPPFLAGS="-DHPL_CALL_CBLAS -DHPL_DETAILED_TIMING"
make -j
```

## 运行过程

```bash
cd testing/ptest
mpiexec -n 4 ../xhpl_ai
```
