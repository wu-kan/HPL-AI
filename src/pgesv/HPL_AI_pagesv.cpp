/*
 * Include files
 */
#include "hpl_ai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPL_AI_pagesv
(
   HPL_T_grid *                     GRID,
   HPL_T_palg *                     ALGO,
   HPL_AI_T_pmat *                  A
)
#else
void HPL_AI_pdgesv
( GRID, ALGO, A )
   HPL_T_grid *                     GRID;
   HPL_T_palg *                     ALGO;
   HPL_AI_T_pmat *                  A;
#endif
{
/* ..
 * .. Executable Statements ..
 */
   if( A->n <= 0 ) return;

   A->info = 0;

   if( ( ALGO->depth == 0 ) || ( GRID->npcol == 1 ) )
   {
      HPL_AI_pagesv0(  GRID, ALGO, A );
   }
   else
   {
      HPL_AI_pagesvK2( GRID, ALGO, A );
   }
/*
 * Solve upper triangular system
 */
   if( A->info == 0 ) HPL_AI_patrsv( GRID, A );
/*
 * End of HPL_AI_pagesv
 */
}

#ifdef __cplusplus
}
#endif
