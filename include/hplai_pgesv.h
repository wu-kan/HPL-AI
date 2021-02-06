#ifndef HPLAI_PGESV_H
#define HPLAI_PGESV_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include <math.h>
#include "hplai.h"

typedef struct HPLAI_S_pmat
{
#ifdef HPL_CALL_VSIPL
   vsip_block_d        * block;
#endif
   HPLAI_T_AFLOAT     * A;            /* pointer to local piece of A */
   HPLAI_T_AFLOAT     * X;             /* pointer to solution vector */
   int                 n;                      /* global problem size */
   int                 nb;                         /* blocking factor */
   int                 ld;                 /* local leading dimension */
   int                 mp;                    /* local number of rows */
   int                 nq;                 /* local number of columns */
   int                 info;                    /* computational flag */
} HPLAI_T_pmat;

void                             HPLAI_pdgesv
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_T_pmat *
) );

void                             HPLAI_pagesv0
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPLAI_T_pmat *
) );
void                             HPLAI_pagesvK1
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPLAI_T_pmat *
) );
void                             HPLAI_pagesvK2
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPLAI_T_pmat *
) );
void                             HPLAI_pagesv
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPLAI_T_pmat *
) );
 
void                             HPLAI_patrsv
STDC_ARGS( (
   HPL_T_grid *,
   HPLAI_T_pmat *
) );

#endif
/*
 * End of hplai_pgesv.h
 */
