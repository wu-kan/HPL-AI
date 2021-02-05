#ifndef HPL_AI_PGESV_H
#define HPL_AI_PGESV_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include <math.h>
#include "hpl_ai.h"

typedef struct HPL_AI_S_pmat
{
#ifdef HPL_CALL_VSIPL
   vsip_block_d        * block;
#endif
   HPL_AI_T_afloat     * A;            /* pointer to local piece of A */
   HPL_AI_T_afloat     * X;             /* pointer to solution vector */
   int                 n;                      /* global problem size */
   int                 nb;                         /* blocking factor */
   int                 ld;                 /* local leading dimension */
   int                 mp;                    /* local number of rows */
   int                 nq;                 /* local number of columns */
   int                 info;                    /* computational flag */
} HPL_AI_T_pmat;

void                             HPL_AI_pdgesv
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_T_pmat *
) );

void                             HPL_AI_pagesv0
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_AI_T_pmat *
) );
void                             HPL_AI_pagesvK1
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_AI_T_pmat *
) );
void                             HPL_AI_pagesvK2
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_AI_T_pmat *
) );
void                             HPL_AI_pagesv
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_AI_T_pmat *
) );
 
void                             HPL_AI_patrsv
STDC_ARGS( (
   HPL_T_grid *,
   HPL_AI_T_pmat *
) );

#endif
/*
 * End of hpl_ai_pgesv.h
 */
