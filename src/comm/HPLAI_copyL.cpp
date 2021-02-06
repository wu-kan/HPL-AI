/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPLAI_copyL
(
   HPLAI_T_panel *                    PANEL
)
#else
void HPLAI_copyL
( PANEL )
   HPLAI_T_panel *                    PANEL;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_copyL copies  the  panel of columns, the L1 replicated submatrix,
 * the pivot array  and  the info scalar into a contiguous workspace for
 * later broadcast.
 *  
 * The copy of this panel  into  a contiguous buffer  can be enforced by
 * specifying -DHPL_COPY_L in the architecture specific Makefile.
 *
 * Arguments
 * =========
 *
 * PANEL   (input/output)                HPLAI_T_panel *
 *         On entry,  PANEL  points to the  current panel data structure
 *         being broadcast.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   int                        jb, lda;
/* ..
 * .. Executable Statements ..
 */
   if( PANEL->grid->mycol == PANEL->pcol )
   {
      jb = PANEL->jb; lda = PANEL->lda;
 
      if( PANEL->grid->myrow == PANEL->prow )
      {
         HPLAI_alacpy( PANEL->mp-jb, jb, Mptr( PANEL->A, jb, -jb, lda ),
                     lda, PANEL->L2, PANEL->ldl2 );
      }
      else
      {
         HPLAI_alacpy( PANEL->mp,    jb, Mptr( PANEL->A,  0, -jb, lda ),
                     lda, PANEL->L2, PANEL->ldl2 );
      }
   }
/*
 * End of HPLAI_copyL
 */
}

#ifdef __cplusplus
}
#endif
