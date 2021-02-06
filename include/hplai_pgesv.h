#ifndef HPLAI_PGESV_H
#define HPLAI_PGESV_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_pgesv.h"

#include "hplai_misc.h"
#include "hplai_blas.h"
#include "hplai_auxil.h"

#include "hplai_grid.h"
#include "hplai_comm.h"
#include "hplai_panel.h"
#include "hplai_pfact.h"

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
