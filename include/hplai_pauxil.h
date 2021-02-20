#ifndef HPLAI_PAUXIL_H
#define HPLAI_PAUXIL_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_pauxil.h"

#include "hplai_misc.h"
#include "hplai_auxil.h"

#include "hplai_grid.h"

#define HPLAI_indxg2lp HPL_indxg2lp
#define HPLAI_indxg2l HPL_indxg2l
#define HPLAI_indxg2p HPL_indxg2p
#define HPLAI_indxl2g HPL_indxl2g
#define HPLAI_infog2l HPL_infog2l
#define HPLAI_numroc HPL_numroc
#define HPLAI_numrocI HPL_numrocI

void                             HPLAI_alaswp00N
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const int *
) );
void                             HPLAI_alaswp10N
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const int *
) );
void                             HPLAI_alaswp01N
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const int *,
   const int *
) );
void                             HPLAI_alaswp01T
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const int *,
   const int *
) );
void                             HPLAI_alaswp02N
STDC_ARGS( (
   const int,
   const int,
   const HPLAI_T_AFLOAT *,
   const int,
   HPLAI_T_AFLOAT *,
   HPLAI_T_AFLOAT *,
   const int,
   const int *,
   const int *
) );
void                             HPLAI_alaswp03N
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const HPLAI_T_AFLOAT *,
   const HPLAI_T_AFLOAT *,
   const int
) );
void                             HPLAI_alaswp03T
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const HPLAI_T_AFLOAT *,
   const HPLAI_T_AFLOAT *,
   const int
) );
void                             HPLAI_alaswp04N
STDC_ARGS( (
   const int,
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const HPLAI_T_AFLOAT *,
   const HPLAI_T_AFLOAT *,
   const int,
   const int *,
   const int *
) );
void                             HPLAI_alaswp04T
STDC_ARGS( (
   const int,
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const HPLAI_T_AFLOAT *,
   const HPLAI_T_AFLOAT *,
   const int,
   const int *,
   const int *
) );
void                             HPLAI_alaswp05N
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const HPLAI_T_AFLOAT *,
   const int,
   const int *,
   const int *
) );
void                             HPLAI_alaswp05T
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const HPLAI_T_AFLOAT *,
   const int,
   const int *,
   const int *
) );
void                             HPLAI_alaswp06N
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const int *
) );
void                             HPLAI_alaswp06T
STDC_ARGS( (
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const int *
) );

#define HPLAI_pabort HPL_pabort
#define HPLAI_pwarn HPL_pwarn
void                             HPLAI_palaprnt
STDC_ARGS( (
   const HPLAI_T_grid *,
   const int,
   const int,
   const int,
   HPLAI_T_AFLOAT *,
   const int,
   const int,
   const int,
   const char *
) );
HPLAI_T_AFLOAT                           HPLAI_palamch
STDC_ARGS( (
   MPI_Comm,
   const HPLAI_T_MACH
) );
HPLAI_T_AFLOAT                           HPLAI_palange
STDC_ARGS( (
   const HPLAI_T_grid *,
   const HPLAI_T_NORM,
   const int,
   const int,
   const int,
   const HPLAI_T_AFLOAT *,
   const int
) );

#endif
/*
 * End of hpl_pauxil.h
 */
