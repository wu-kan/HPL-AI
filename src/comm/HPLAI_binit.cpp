/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
   int HPLAI_binit(
       HPLAI_T_panel *PANEL)
#else
int HPLAI_binit(PANEL)
    HPLAI_T_panel *PANEL;
#endif
   {
      /*
 * .. Local Variables ..
 */
      int ierr;
      HPL_T_TOP top;
      /* ..
 * .. Executable Statements ..
 */
      if (PANEL->grid->npcol <= 1)
         return (HPL_SUCCESS);
      /*
 * Retrieve the selected virtual broadcast topology
 */
      top = PANEL->algo->btopo;

      switch (top)
      {
      case HPL_1RING_M:
         ierr = HPLAI_binit_1rinM(PANEL);
         break;
      case HPL_1RING:
         ierr = HPLAI_binit_1ring(PANEL);
         break;
      case HPL_2RING_M:
         ierr = HPLAI_binit_2rinM(PANEL);
         break;
      case HPL_2RING:
         ierr = HPLAI_binit_2ring(PANEL);
         break;
      case HPL_BLONG_M:
         ierr = HPLAI_binit_blonM(PANEL);
         break;
      case HPL_BLONG:
         ierr = HPLAI_binit_blong(PANEL);
         break;
      default:
         ierr = HPL_SUCCESS;
      }

      return (ierr);
      /*
 * End of HPLAI_binit
 */
   }

#ifdef __cplusplus
}
#endif
