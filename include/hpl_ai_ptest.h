#ifndef HPL_AI_PTEST_H
#define HPL_AI_PTEST_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_ptest.h"
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_AI_pdinfo
STDC_ARGS( (
   HPL_T_test *,
   int *,
   int *,
   int *,
   int *,
   HPL_T_ORDER *,
   int *,
   int *,
   int *,
   int *,
   HPL_T_FACT *,
   int *,
   int *,
   int *,
   int *,
   int *,
   HPL_T_FACT *,
   int *,
   HPL_T_TOP *,
   int *,
   int *,
   HPL_T_SWAP *,
   int *,
   int *,
   int *,
   int *,
   int *
) );
void                             HPL_AI_pdtest
STDC_ARGS( (
   HPL_T_test *,
   HPL_T_grid *,
   HPL_T_palg *,
   const int,
   const int
) );

#endif
/*
 * End of hpl_ai_ptest.h
 */
