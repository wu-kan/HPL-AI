/*
 * Include files
 */
#include "hplai.hh"

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
int HPLAI_recv
(
   HPLAI_T_AFLOAT *                         RBUF,
   int                              RCOUNT,
   int                              SRC,
   int                              RTAG,
   MPI_Comm                         COMM
)
#else
int HPLAI_recv
( RBUF, RCOUNT, SRC, RTAG, COMM )
   HPLAI_T_AFLOAT *                         RBUF;
   int                              RCOUNT;
   int                              SRC;
   int                              RTAG;
   MPI_Comm                         COMM;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_recv is a simple wrapper around  MPI_Recv.  Its  main  purpose is
 * to  allow for some  experimentation / tuning  of this simple routine.
 * Successful  completion  is  indicated  by  the  returned  error  code
 * HPL_SUCCESS.  In the case of messages of length less than or equal to
 * zero, this function returns immediately.
 *
 * Arguments
 * =========
 *
 * RBUF    (local output)                HPLAI_T_AFLOAT *
 *         On entry, RBUF specifies the starting address of buffer to be
 *         received.
 *
 * RCOUNT  (local input)                 int
 *         On entry,  RCOUNT  specifies  the number  of HPLAI_T_AFLOAT precision
 *         entries in RBUF. RCOUNT must be at least zero.
 *
 * SRC     (local input)                 int
 *         On entry, SRC  specifies the rank of the  sending  process in
 *         the communication space defined by COMM.
 *
 * RTAG    (local input)                 int
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
   MPI_Status                 status;
#ifdef HPL_USE_MPI_DATATYPE
   MPI_Datatype               type;
#endif
   int                        ierr;
/* ..
 * .. Executable Statements ..
 */
   if( RCOUNT <= 0 ) return( HPL_SUCCESS );

#ifdef HPL_USE_MPI_DATATYPE
   ierr =      MPI_Type_contiguous( RCOUNT, HPLAI_MPI_AFLOAT, &type );
   if( ierr == MPI_SUCCESS )
      ierr =   MPI_Type_commit( &type );
   if( ierr == MPI_SUCCESS )
      ierr =   MPI_Recv( (void *)(RBUF), 1, type, SRC, RTAG, COMM,
                         &status );
   if( ierr == MPI_SUCCESS )
      ierr =   MPI_Type_free( &type );
#else
   ierr = MPI_Recv( (void *)(RBUF), RCOUNT, HPLAI_MPI_AFLOAT, SRC, RTAG,
                    COMM, &status );
#endif
   return( ( ierr == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE ) );
/*
 * End of HPLAI_recv
 */
}

#ifdef __cplusplus
}
#endif
