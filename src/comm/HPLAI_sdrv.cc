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
int HPLAI_sdrv
(
   HPLAI_T_AFLOAT *                         SBUF,
   int                              SCOUNT,
   int                              STAG,
   HPLAI_T_AFLOAT *                         RBUF,
   int                              RCOUNT,
   int                              RTAG,
   int                              PARTNER,
   MPI_Comm                         COMM
)
#else
int HPLAI_sdrv
( SBUF, SCOUNT, STAG, RBUF, RCOUNT, RTAG, PARTNER, COMM )
   HPLAI_T_AFLOAT *                         SBUF;
   int                              SCOUNT;
   int                              STAG;
   HPLAI_T_AFLOAT *                         RBUF;
   int                              RCOUNT;
   int                              RTAG;
   int                              PARTNER;
   MPI_Comm                         COMM;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_sdrv is a simple wrapper around MPI_Sendrecv. Its main purpose is
 * to allow for some experimentation and tuning of this simple function.
 * Messages  of  length  less than  or  equal to zero  are not sent  nor
 * received.  Successful completion  is  indicated by the returned error
 * code HPL_SUCCESS.
 *
 * Arguments
 * =========
 *
 * SBUF    (local input)                 HPLAI_T_AFLOAT *
 *         On entry, SBUF specifies the starting address of buffer to be
 *         sent.
 *
 * SCOUNT  (local input)                 int
 *         On entry,  SCOUNT  specifies  the number  of HPLAI_T_AFLOAT precision
 *         entries in SBUF. SCOUNT must be at least zero.
 *
 * STAG    (local input)                 int
 *         On entry,  STAG  specifies the message tag to be used for the
 *         sending communication operation.
 *
 * RBUF    (local output)                HPLAI_T_AFLOAT *
 *         On entry, RBUF specifies the starting address of buffer to be
 *         received.
 *
 * RCOUNT  (local input)                 int
 *         On entry,  RCOUNT  specifies  the number  of HPLAI_T_AFLOAT precision
 *         entries in RBUF. RCOUNT must be at least zero.
 *
 * RTAG    (local input)                 int
 *         On entry,  RTAG  specifies the message tag to be used for the
 *         receiving communication operation.
 *
 * PARTNER (local input)                 int
 *         On entry,  PARTNER  specifies  the rank of the  collaborative
 *         process in the communication space defined by COMM.
 *
 * COMM    (local input)                 MPI_Comm
 *         The MPI communicator identifying the communication space.
 *
 * ---------------------------------------------------------------------
 */
   if(SCOUNT && RCOUNT)
      return( ( MPI_Sendrecv((const void *)(SBUF), SCOUNT, HPLAI_MPI_AFLOAT, PARTNER, STAG, (void *)(RBUF), RCOUNT, HPLAI_MPI_AFLOAT, PARTNER, RTAG, COMM, MPI_STATUS_IGNORE) == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE ) );
/*
 * .. Local Variables ..
 */
#ifdef HPL_USE_MPI_DATATYPE
   MPI_Datatype               type[2];
#endif
   MPI_Request                request;
   MPI_Status                 status;
   int                        ierr;
/* ..
 * .. Executable Statements ..
 */
   if( RCOUNT > 0 )
   {
      if( SCOUNT > 0 )
      {
#ifdef HPL_USE_MPI_DATATYPE
/*
 * Post asynchronous receive
 */
         ierr =      MPI_Type_contiguous( RCOUNT, HPLAI_MPI_AFLOAT, &type[0] );
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Type_commit( &type[0] );
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Irecv( (void *)(RBUF), 1, type[0], PARTNER,
                                RTAG, COMM, &request );
/*
 * Blocking send
 */
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Type_contiguous( SCOUNT, HPLAI_MPI_AFLOAT, &type[1] );
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Type_commit( &type[1] );
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Send( (void *)(SBUF), 1, type[1], PARTNER,
                               STAG, COMM );
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Type_free( &type[1] );
/*
 * Wait for the receive to complete
 */
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Wait( &request, &status );
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Type_free( &type[0] );
#else
/*
 * Post asynchronous receive
 */
         ierr =      MPI_Irecv( (void *)(RBUF), RCOUNT, HPLAI_MPI_AFLOAT,
                                PARTNER, RTAG, COMM, &request );
/*
 * Blocking send
 */
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Send( (void *)(SBUF), SCOUNT, HPLAI_MPI_AFLOAT,
                               PARTNER, STAG, COMM );
/*
 * Wait for the receive to complete
 */
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Wait( &request, &status );
#endif
      }
      else
      {
/*
 * Blocking receive
 */
#ifdef HPL_USE_MPI_DATATYPE
         ierr =      MPI_Type_contiguous( RCOUNT, HPLAI_MPI_AFLOAT, &type[0] );
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Type_commit( &type[0] );
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Recv( (void *)(RBUF), 1, type[0], PARTNER, RTAG,
                               COMM, &status );
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Type_free( &type[0] );
#else
         ierr =      MPI_Recv( (void *)(RBUF), RCOUNT, HPLAI_MPI_AFLOAT,
                               PARTNER, RTAG, COMM, &status );
#endif
      }
   }
   else if( SCOUNT > 0 )
   {
/*
 * Blocking send
 */
#ifdef HPL_USE_MPI_DATATYPE
      ierr =      MPI_Type_contiguous( SCOUNT, HPLAI_MPI_AFLOAT, &type[1] );
      if( ierr == MPI_SUCCESS )
         ierr =   MPI_Type_commit( &type[1] );
      if( ierr == MPI_SUCCESS )
         ierr =   MPI_Send( (void *)(SBUF), 1, type[1], PARTNER, STAG,
                          COMM );
      if( ierr == MPI_SUCCESS )
         ierr =   MPI_Type_free( &type[1] ) );
#else
      ierr =      MPI_Send( (void *)(SBUF), SCOUNT, HPLAI_MPI_AFLOAT, PARTNER,
                            STAG, COMM );
#endif
   }
   else { ierr = MPI_SUCCESS; }

   return( ( ierr == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE ) );
/*
 * End of HPLAI_sdrv
 */
}

#ifdef __cplusplus
}
#endif
