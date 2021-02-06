#ifndef HPLAI_COMM_H
#define HPLAI_COMM_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hplai.h"

/*
 * ---------------------------------------------------------------------
 * comm function prototypes
 * ---------------------------------------------------------------------
 */
int                              HPLAI_send
STDC_ARGS( (
   HPLAI_T_AFLOAT *,
   int,
   int,
   int,
   MPI_Comm
) );
int                              HPLAI_recv
STDC_ARGS( (
   HPLAI_T_AFLOAT *,
   int,
   int,
   int,
   MPI_Comm
) );

#endif
/*
 * End of hplai_comm.h
 */
