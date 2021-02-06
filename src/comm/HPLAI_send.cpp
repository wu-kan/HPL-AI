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
   int HPLAI_send(
       HPLAI_T_AFLOAT *SBUF,
       int SCOUNT,
       int DEST,
       int STAG,
       MPI_Comm COMM)
#else
int HPLAI_send(SBUF, SCOUNT, DEST, STAG, COMM)
    HPLAI_T_AFLOAT *SBUF;
int SCOUNT;
int DEST;
int STAG;
MPI_Comm COMM;
#endif
   {
/*
 * .. Local Variables ..
 */
#ifdef HPL_USE_MPI_DATATYPE
      MPI_Datatype type;
#endif
      int ierr;
      /* ..
 * .. Executable Statements ..
 */
      if (SCOUNT <= 0)
         return (HPL_SUCCESS);

#ifdef HPL_USE_MPI_DATATYPE
      ierr = MPI_Type_contiguous(1LL * sizeof(HPLAI_T_AFLOAT) * SCOUNT, MPI_BYTE, &type);
      if (ierr == MPI_SUCCESS)
         ierr = MPI_Type_commit(&type);
      if (ierr == MPI_SUCCESS)
         ierr = MPI_Send((void *)(SBUF), 1, type, DEST, STAG, COMM);
      if (ierr == MPI_SUCCESS)
         ierr = MPI_Type_free(&type);
#else
   ierr = MPI_Send((void *)(SBUF), 1LL * sizeof(HPLAI_T_AFLOAT) * SCOUNT, MPI_BYTE, DEST, STAG, COMM);
#endif
      return ((ierr == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE));
      /*
 * End of HPL_send
 */
   }

#ifdef __cplusplus
}
#endif
