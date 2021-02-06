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
