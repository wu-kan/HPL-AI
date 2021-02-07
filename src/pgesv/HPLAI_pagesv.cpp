/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPLAI_pagesv
(
   HPLAI_T_grid *                     GRID,
   HPLAI_T_palg *                   ALGO,
   HPLAI_T_pmat *                   A
)
#else
void HPLAI_pdgesv
( GRID, ALGO, A )
   HPLAI_T_grid *                     GRID;
   HPLAI_T_palg *                   ALGO;
   HPLAI_T_pmat *                   A;
#endif
{
/* ..
 * .. Executable Statements ..
 */
   if( A->n <= 0 ) return;

   A->info = 0;

   if( ( ALGO->depth == 0 ) || ( GRID->npcol == 1 ) )
   {
      HPLAI_pagesv0(  GRID, ALGO, A );
   }
   else
   {
      HPLAI_pagesvK2( GRID, ALGO, A );
   }
/*
 * Solve upper triangular system
 */
   if( A->info == 0 ) HPLAI_patrsv( GRID, A );
/*
 * End of HPLAI_pagesv
 */
}

#ifdef __cplusplus
}
#endif
