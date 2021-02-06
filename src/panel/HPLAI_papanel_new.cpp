/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPLAI_papanel_new
(
   HPL_T_grid *                     GRID,
   HPL_T_palg *                     ALGO,
   const int                        M,
   const int                        N,
   const int                        JB,
   HPLAI_T_pmat *                  A,
   const int                        IA,
   const int                        JA,
   const int                        TAG,
   HPLAI_T_panel * *               PANEL
)
#else
void HPLAI_papanel_new
( GRID, ALGO, M, N, JB, A, IA, JA, TAG, PANEL )
   HPL_T_grid *                     GRID;
   HPL_T_palg *                     ALGO;
   const int                        M;
   const int                        N;
   const int                        JB;
   HPLAI_T_pmat *                  A;
   const int                        IA;
   const int                        JA;
   const int                        TAG;
   HPLAI_T_panel * *               PANEL;
#endif
{
/*
 * .. Local Variables ..
 */
   HPLAI_T_panel                * p = NULL;
/* ..
 * .. Executable Statements ..
 */
/*
 * Allocate the panel structure - Check for enough memory
 */
   if( !( p = (HPLAI_T_panel *)malloc( sizeof( HPLAI_T_panel ) ) ) )
   {
      HPL_pabort( __LINE__, "HPLAI_papanel_new", "Memory allocation failed" );
   }

   HPLAI_papanel_init( GRID, ALGO, M, N, JB, A, IA, JA, TAG, p );
   *PANEL = p;
/*
 * End of HPLAI_papanel_new
 */
}

#ifdef __cplusplus
}
#endif
