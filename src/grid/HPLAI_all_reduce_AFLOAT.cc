/*
 * Include files
 */
#include "hplai.hh"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
int HPLAI_all_reduce_AFLOAT
(
   void *                           BUFFER,
   const int                        COUNT,
   const HPLAI_T_OP_AFLOAT          OP,
   MPI_Comm                         COMM
)
#else
int HPLAI_all_reduce_AFLOAT
( BUFFER, COUNT, DTYPE, OP, COMM )
   void *                           BUFFER;
   const int                        COUNT;
   const HPLAI_T_OP_AFLOAT          OP;
   MPI_Comm                         COMM;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_all_reduce_AFLOAT performs   a   global   reduce  operation  across  all
 * processes of a group leaving the results on all processes.
 *
 * Arguments
 * =========
 *
 * BUFFER  (local input/global output)   void *
 *         On entry,  BUFFER  points to  the  buffer to be combined.  On
 *         exit, this array contains the combined data and  is identical
 *         on all processes in the group.
 *
 * COUNT   (global input)                const int
 *         On entry,  COUNT  indicates the number of entries in  BUFFER.
 *         COUNT must be at least zero.
 *
 * OP      (global input)                const HPLAI_T_OP_AFLOAT 
 *         On entry, OP is a pointer to the local combine function.
 *
 * COMM    (global/local input)          MPI_Comm
 *         The MPI communicator identifying the process collection.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   int                        hplerr;
/* ..
 * .. Executable Statements ..
 */
   hplerr = HPLAI_reduce_AFLOAT(   BUFFER, COUNT, OP, 0, COMM );
   if( hplerr != MPI_SUCCESS ) return( hplerr );
   return( HPLAI_broadcast_AFLOAT( BUFFER, COUNT,     0, COMM ) );
/*
 * End of HPLAI_all_reduce_AFLOAT
 */
}

#ifdef __cplusplus
}
#endif
