/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
int HPLAI_papanel_disp
(
   HPLAI_T_panel * *                  PANEL
)
#else
int HPLAI_papanel_disp
( PANEL )
   HPLAI_T_panel * *                  PANEL;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_papanel_disp deallocates  the  panel  structure  and  resources  and
 * stores the error code returned by the panel factorization.
 *
 * Arguments
 * =========
 *
 * PANEL   (local input/output)          HPLAI_T_panel * *
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
   mpierr = HPLAI_papanel_free( *PANEL );
   if( *PANEL ) free( *PANEL );
   *PANEL = NULL;

   return( mpierr );
/*
 * End of HPLAI_papanel_disp
 */
}

#ifdef __cplusplus
}
#endif
