/*
 * Include files
 */
#include "hplai.hh"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
   int HPLAI_bcast(
       HPLAI_T_panel *PANEL,
       int *IFLAG)
#else
int HPLAI_bcast(PANEL, IFLAG)
    HPLAI_T_panel *PANEL;
int *IFLAG;
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
      if (PANEL == NULL)
      {
         *IFLAG = HPL_SUCCESS;
         return (HPL_SUCCESS);
      }
      if (PANEL->grid->npcol <= 1)
      {
         *IFLAG = HPL_SUCCESS;
         return (HPL_SUCCESS);
      }
      /*
 * Retrieve the selected virtual broadcast topology
 */
      top = PANEL->algo->btopo;

      switch (top)
      {
      case HPL_1RING_M:
         ierr = HPLAI_bcast_1rinM(PANEL, IFLAG);
         break;
      case HPL_1RING:
         ierr = HPLAI_bcast_1ring(PANEL, IFLAG);
         break;
      case HPL_2RING_M:
         ierr = HPLAI_bcast_2rinM(PANEL, IFLAG);
         break;
      case HPL_2RING:
         ierr = HPLAI_bcast_2ring(PANEL, IFLAG);
         break;
      case HPL_BLONG_M:
         ierr = HPLAI_bcast_blonM(PANEL, IFLAG);
         break;
      case HPL_BLONG:
         ierr = HPLAI_bcast_blong(PANEL, IFLAG);
         break;
      default:
         ierr = HPL_SUCCESS;
      }

      return (ierr);
      /*
 * End of HPLAI_bcast
 */
   }

#ifdef __cplusplus
}
#endif
