
AM_CPPFLAGS = -I$(top_srcdir)/include

xhpl_LDADD = ../src/libhpl.a

xhpl_ai_LDADD = ../src/libhpl_ai.a ../src/libhpl.a

bin_PROGRAMS = xhpl_ai xhpl

xhpl_ai_SOURCES =  \
matgen/HPL_jumpit.c matgen/HPL_rand.c matgen/HPL_setran.c matgen/HPL_xjumpm.c \
matgen/HPL_lmul.c matgen/HPL_ladd.c \
pmatgen/HPLAI_pdmatgen.cc \
ptest/HPLAI_pddriver.cc ptest/HPLAI_pdinfo.cc ptest/HPLAI_pdtest.cc \
ptimer/HPL_ptimer.c ptimer/HPL_ptimer_cputime.c ptimer/HPL_ptimer_walltime.c

xhpl_SOURCES =  \
matgen/HPL_jumpit.c matgen/HPL_rand.c matgen/HPL_setran.c matgen/HPL_xjumpm.c \
matgen/HPL_lmul.c matgen/HPL_ladd.c \
pmatgen/HPL_pdmatgen.c \
ptest/HPL_pddriver.c ptest/HPL_pdinfo.c ptest/HPL_pdtest.c \
ptimer/HPL_ptimer.c ptimer/HPL_ptimer_cputime.c ptimer/HPL_ptimer_walltime.c
