/*
 * Include files
 */
#include "hpl_ai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPL_AI_papanel_new
(
   HPL_T_grid *                     GRID,
   HPL_T_palg *                     ALGO,
   const int                        M,
   const int                        N,
   const int                        JB,
   HPL_AI_T_pmat *                  A,
   const int                        IA,
   const int                        JA,
   const int                        TAG,
   HPL_AI_T_panel * *               PANEL
)
#else
void HPL_AI_papanel_new
( GRID, ALGO, M, N, JB, A, IA, JA, TAG, PANEL )
   HPL_T_grid *                     GRID;
   HPL_T_palg *                     ALGO;
   const int                        M;
   const int                        N;
   const int                        JB;
   HPL_AI_T_pmat *                  A;
   const int                        IA;
   const int                        JA;
   const int                        TAG;
   HPL_AI_T_panel * *               PANEL;
#endif
{
/*
 * .. Local Variables ..
 */
   HPL_AI_T_panel                * p = NULL;
/* ..
 * .. Executable Statements ..
 */
/*
 * Allocate the panel structure - Check for enough memory
 */
   if( !( p = (HPL_AI_T_panel *)malloc( sizeof( HPL_AI_T_panel ) ) ) )
   {
      HPL_pabort( __LINE__, "HPL_AI_papanel_new", "Memory allocation failed" );
   }

   HPL_AI_papanel_init( GRID, ALGO, M, N, JB, A, IA, JA, TAG, p );
   *PANEL = p;
/*
 * End of HPL_AI_papanel_new
 */
}

#ifdef __cplusplus
}
#endif
