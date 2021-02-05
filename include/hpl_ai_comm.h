#ifndef HPL_AI_COMM_H
#define HPL_AI_COMM_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_ai.h"

/*
 * ---------------------------------------------------------------------
 * comm function prototypes
 * ---------------------------------------------------------------------
 */
int                              HPL_AI_send
STDC_ARGS( (
   HPL_AI_T_afloat *,
   int,
   int,
   int,
   MPI_Comm
) );
int                              HPL_AI_recv
STDC_ARGS( (
   HPL_AI_T_afloat *,
   int,
   int,
   int,
   MPI_Comm
) );

#endif
/*
 * End of hpl_ai_comm.h
 */
