#ifndef HPLAI_COMM_H
#define HPLAI_COMM_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_comm.h"
#include "hplai_misc.h"
#include "hplai_panel.h"

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
int                              HPLAI_sdrv
STDC_ARGS( (
   HPLAI_T_AFLOAT *,
   int,
   int,
   HPLAI_T_AFLOAT *,
   int,
   int,
   int,
   MPI_Comm
) );
int                              HPLAI_binit
STDC_ARGS( (
   HPLAI_T_panel *
) );
int                              HPLAI_bcast
STDC_ARGS( (
   HPLAI_T_panel *,
   int *
) );
int                              HPLAI_bwait
STDC_ARGS( (
   HPLAI_T_panel *
) );
int                              HPLAI_packL
STDC_ARGS( (
   HPLAI_T_panel *,
   const int,
   const int,
   const int
) );
 
int HPLAI_binit_1ring STDC_ARGS( ( HPLAI_T_panel *        ) );
int HPLAI_bcast_1ring STDC_ARGS( ( HPLAI_T_panel *, int * ) );
int HPLAI_bwait_1ring STDC_ARGS( ( HPLAI_T_panel *        ) );
 
int HPLAI_binit_1rinM STDC_ARGS( ( HPLAI_T_panel *        ) );
int HPLAI_bcast_1rinM STDC_ARGS( ( HPLAI_T_panel *, int * ) );
int HPLAI_bwait_1rinM STDC_ARGS( ( HPLAI_T_panel *        ) );
 
int HPLAI_binit_2ring STDC_ARGS( ( HPLAI_T_panel *        ) );
int HPLAI_bcast_2ring STDC_ARGS( ( HPLAI_T_panel *, int * ) );
int HPLAI_bwait_2ring STDC_ARGS( ( HPLAI_T_panel *        ) );
 
int HPLAI_binit_2rinM STDC_ARGS( ( HPLAI_T_panel *        ) );
int HPLAI_bcast_2rinM STDC_ARGS( ( HPLAI_T_panel *, int * ) );
int HPLAI_bwait_2rinM STDC_ARGS( ( HPLAI_T_panel *        ) );
 
int HPLAI_binit_blong STDC_ARGS( ( HPLAI_T_panel *        ) );
int HPLAI_bcast_blong STDC_ARGS( ( HPLAI_T_panel *, int * ) );
int HPLAI_bwait_blong STDC_ARGS( ( HPLAI_T_panel *        ) );
 
int HPLAI_binit_blonM STDC_ARGS( ( HPLAI_T_panel *        ) );
int HPLAI_bcast_blonM STDC_ARGS( ( HPLAI_T_panel *, int * ) );
int HPLAI_bwait_blonM STDC_ARGS( ( HPLAI_T_panel *        ) );

#endif
/*
 * End of hplai_comm.h
 */
