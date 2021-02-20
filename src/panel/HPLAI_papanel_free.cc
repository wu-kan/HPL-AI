#include "hplai.hh"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
int HPLAI_papanel_free
(
   HPLAI_T_panel *                    PANEL
)
#else
int HPLAI_papanel_free
( PANEL )
   HPLAI_T_panel *                    PANEL;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_pdpanel_free deallocates  the panel resources  and  stores the error
 * code returned by the panel factorization.
 *
 * Arguments
 * =========
 *
 * PANEL   (local input/output)          HPL_T_panel *
 *         On entry,  PANEL  points  to  the  panel data  structure from
 *         which the resources should be deallocated.
 *
 * ---------------------------------------------------------------------
 */ 
/* ..
 * .. Executable Statements ..
 */
   if( PANEL->pmat->info == 0 ) PANEL->pmat->info = *(PANEL->DINFO);
#ifdef HPL_CALL_VSIPL
/*
 * Release the blocks
 */
   (void) vsip_blockrelease_d( PANEL->L1block, VSIP_TRUE );
   (void) vsip_blockrelease_d( PANEL->L2block, VSIP_TRUE );
   if( PANEL->grid->nprow > 1 )
      (void) vsip_blockrelease_d( PANEL->Ublock,  VSIP_TRUE );
/*
 * Destroy blocks
 */
   vsip_blockdestroy_d( PANEL->L1block );
   vsip_blockdestroy_d( PANEL->L2block );
   if( PANEL->grid->nprow > 1 )
      vsip_blockdestroy_d( PANEL->Ublock );
#endif

   if( PANEL->WORK  ) free( PANEL->WORK  );
   if( PANEL->IWORK ) free( PANEL->IWORK );

   return( MPI_SUCCESS );
/*
 * End of HPLAI_papanel_free
 */
}

#ifdef __cplusplus
}
#endif
