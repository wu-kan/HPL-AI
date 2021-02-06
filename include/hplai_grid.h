#ifndef HPLAI_GRID_H
#define HPLAI_GRID_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hplai.h"

int                              HPLAI_broadcast_afloat
STDC_ARGS( (
   void *,
   const int,
   const int,
   MPI_Comm
) );

#endif
/*
 * End of hplai_grid.h
 */
