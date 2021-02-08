#ifndef HPLAI_AUXIL_H
#define HPLAI_AUXIL_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_auxil.h"
#include "hplai_misc.h"
#include "hplai_blas.h"

#define HPLAI_NORM_A HPL_NORM_A
#define HPLAI_NORM_1 HPL_NORM_1
#define HPLAI_NORM_I HPL_NORM_I
#define HPLAI_T_NORM HPL_T_NORM

#define HPLAI_MACH_EPS HPL_MACH_EPS
#define HPLAI_MACH_SFMIN HPL_MACH_SFMIN
#define HPLAI_MACH_BASE HPL_MACH_BASE
#define HPLAI_MACH_PREC HPL_MACH_PREC
#define HPLAI_MACH_MLEN HPL_MACH_MLEN
#define HPLAI_MACH_RND HPL_MACH_RND
#define HPLAI_MACH_EMIN HPL_MACH_EMIN
#define HPLAI_MACH_RMIN HPL_MACH_RMIN
#define HPLAI_MACH_EMAX HPL_MACH_EMAX
#define HPLAI_MACH_RMAX HPL_MACH_RMAX
#define HPLAI_T_MACH HPL_T_MACH

#define HPLAI_fprintf HPL_fprintf
#define HPLAI_warn HPL_warn
#define HPLAI_abort HPL_abort

void                             HPLAI_alacpy
STDC_ARGS( (
   const int,
   const int,
   const HPLAI_T_AFLOAT *,
   const int,
   HPLAI_T_AFLOAT *,
   const int
) );
void                             HPLAI_alatcpy
STDC_ARGS( (
   const int,
   const int,
   const HPLAI_T_AFLOAT *,
   const int,
   HPLAI_T_AFLOAT *,
   const int
) );
HPLAI_T_AFLOAT                           HPLAI_alamch
STDC_ARGS( (
   const HPL_T_MACH
) );

#endif
/*
 * End of hpl_auxil.h
 */
