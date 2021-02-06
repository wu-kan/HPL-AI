/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
   int HPLAI_sdrv(
       HPLAI_T_AFLOAT *SBUF,
       int SCOUNT,
       int STAG,
       HPLAI_T_AFLOAT *RBUF,
       int RCOUNT,
       int RTAG,
       int PARTNER,
       MPI_Comm COMM)
#else
int HPLAI_sdrv(SBUF, SCOUNT, STAG, RBUF, RCOUNT, RTAG, PARTNER, COMM) HPLAI_T_AFLOAT *SBUF;
int SCOUNT;
int STAG;
HPLAI_T_AFLOAT *RBUF;
int RCOUNT;
int RTAG;
int PARTNER;
MPI_Comm COMM;
#endif
   {
      return MPI_Sendrecv(
                 (void *)SBUF,
                 1LL * sizeof(HPLAI_T_AFLOAT) * SCOUNT,
                 MPI_BYTE,
                 PARTNER,
                 STAG,
                 (void *)RBUF,
                 1LL * sizeof(HPLAI_T_AFLOAT) * RCOUNT,
                 MPI_BYTE,
                 PARTNER,
                 RTAG,
                 COMM,
                 MPI_STATUS_IGNORE) == MPI_SUCCESS
                 ? HPL_SUCCESS
                 : HPL_FAILURE;
      /*
 * End of HPLAI_sdrv
 */
   }

#ifdef __cplusplus
}
#endif
