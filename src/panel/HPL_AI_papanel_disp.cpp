/*
 * Include files
 */
#include "hpl_ai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
int HPL_AI_papanel_disp
(
   HPL_AI_T_panel * *                  PANEL
)
#else
int HPL_AI_papanel_disp
( PANEL )
   HPL_AI_T_panel * *                  PANEL;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_AI_papanel_disp deallocates  the  panel  structure  and  resources  and
 * stores the error code returned by the panel factorization.
 *
 * Arguments
 * =========
 *
 * PANEL   (local input/output)          HPL_AI_T_panel * *
 *         On entry,  PANEL  points  to  the  address  of the panel data
 *         structure to be deallocated.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   int                        mpierr;
/* ..
 * .. Executable Statements ..
 */
/*
 * Deallocate the panel resources and panel structure
 */
   mpierr = HPL_AI_papanel_free( *PANEL );
   if( *PANEL ) free( *PANEL );
   *PANEL = NULL;

   return( mpierr );
/*
 * End of HPL_AI_papanel_disp
 */
}

#ifdef __cplusplus
}
#endif
