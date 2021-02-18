/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif
/*
 * Do not use  MPI user-defined data types no matter what.  This routine
 * is used for small contiguous messages.
 */
#ifdef HPL_USE_MPI_DATATYPE
#undef HPL_USE_MPI_DATATYPE
#endif

#ifdef STDC_HEADERS
int HPLAI_send
(
   HPLAI_T_AFLOAT *                         SBUF,
   int                              SCOUNT,
   int                              DEST,
   int                              STAG,
   MPI_Comm                         COMM
)
#else
int HPLAI_send
( SBUF, SCOUNT, DEST, STAG, COMM )
   HPLAI_T_AFLOAT *                         SBUF;
   int                              SCOUNT;
   int                              DEST;
   int                              STAG;
   MPI_Comm                         COMM;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_send is a simple wrapper around  MPI_Send.  Its  main  purpose is
 * to  allow for some  experimentation / tuning  of this simple routine.
 * Successful  completion  is  indicated  by  the  returned  error  code
 * MPI_SUCCESS.  In the case of messages of length less than or equal to
 * zero, this function returns immediately.
 *
 * Arguments
 * =========
 *
 * SBUF    (local input)                 HPLAI_T_AFLOAT *
 *         On entry, SBUF specifies the starting address of buffer to be
 *         sent.
 *
 * SCOUNT  (local input)                 int
 *         On entry,  SCOUNT  specifies  the number of  HPLAI_T_AFLOAT precision
 *         entries in SBUF. SCOUNT must be at least zero.
 *
 * DEST    (local input)                 int
 *         On entry, DEST specifies the rank of the receiving process in
 *         the communication space defined by COMM.
 *
 * STAG    (local input)                 int
 *         On entry,  STAG specifies the message tag to be used for this
 *         communication operation.
 *
 * COMM    (local input)                 MPI_Comm
 *         The MPI communicator identifying the communication space.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
#ifdef HPL_USE_MPI_DATATYPE
   MPI_Datatype               type;
#endif
   int                        ierr;
/* ..
 * .. Executable Statements ..
 */
   if( SCOUNT <= 0 ) return( HPL_SUCCESS );

#ifdef HPL_USE_MPI_DATATYPE
   ierr =      MPI_Type_contiguous( SCOUNT, HPLAI_MPI_AFLOAT, &type );
   if( ierr == MPI_SUCCESS )
      ierr =   MPI_Type_commit( &type );
   if( ierr == MPI_SUCCESS )
      ierr =   MPI_Send( (void *)(SBUF), 1, type, DEST, STAG, COMM );
   if( ierr == MPI_SUCCESS )
      ierr =   MPI_Type_free( &type );
#else
   ierr = MPI_Send( (void *)(SBUF), SCOUNT, HPLAI_MPI_AFLOAT, DEST, STAG, COMM );
#endif
   return( ( ierr == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE ) ); 
/*
 * End of HPLAI_send
 */
}

#ifdef __cplusplus
}
#endif
