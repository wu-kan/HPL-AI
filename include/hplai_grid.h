#ifndef HPLAI_GRID_H
#define HPLAI_GRID_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_grid.h"
#include "hplai_misc.h"

/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */
#define HPLAI_ROW_MAJOR HPL_ROW_MAJOR
#define HPLAI_COLUMN_MAJOR HPL_COLUMN_MAJOR
#define HPLAI_T_ORDER HPL_T_ORDER
#define HPLAI_S_grid HPL_S_grid
#define HPLAI_T_grid HPL_T_grid

/*
 * ---------------------------------------------------------------------
 * Data Structures
 * ---------------------------------------------------------------------
 */
typedef void (*HPLAI_T_OP_AFLOAT)
(  const int,       const void *,    void *);

#define HPLAI_grid_init HPL_grid_init
#define HPLAI_grid_exit HPL_grid_exit
#define HPLAI_grid_info HPL_grid_info
#define HPLAI_pnum HPL_pnum
#define HPLAI_barrier HPL_barrier
int                              HPLAI_broadcast_AFLOAT
STDC_ARGS( (
   void *,
   const int,
   const int,
   MPI_Comm
) );
int                              HPLAI_reduce_AFLOAT
STDC_ARGS( (
   void *,
   const int,
   const HPLAI_T_OP_AFLOAT ,
   const int,
   MPI_Comm
) );
int                              HPLAI_all_reduce_AFLOAT
STDC_ARGS( (
   void *,
   const int,
   const HPLAI_T_OP_AFLOAT ,
   MPI_Comm
) );

void                             HPLAI_max_AFLOAT
STDC_ARGS( (
   const int,
   const void *,
   void *
) );
void                             HPLAI_min_AFLOAT
STDC_ARGS( (
   const int,
   const void *,
   void *
) );
void                             HPLAI_sum_AFLOAT
STDC_ARGS( (
   const int,
   const void *,
   void *
) );

#endif
/*
 * End of hplai_grid.h
 */
