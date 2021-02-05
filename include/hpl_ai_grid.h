#ifndef HPL_AI_GRID_H
#define HPL_AI_GRID_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_ai.h"

int                              HPL_AI_broadcast_afloat
STDC_ARGS( (
   void *,
   const int,
   const int,
   MPI_Comm
) );

#endif
/*
 * End of hpl_ai_grid.h
 */
