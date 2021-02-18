/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
int HPLAI_broadcast_AFLOAT(
    void *BUFFER,
    const int COUNT,
    const int ROOT,
    MPI_Comm COMM)
#else
int HPLAI_broadcast_AFLOAT(BUFFER, COUNT, DTYPE, ROOT, COMM) void *BUFFER;
const int COUNT;
const int ROOT;
MPI_Comm COMM;
#endif
{
   return MPI_Bcast(BUFFER, COUNT, HPLAI_MPI_AFLOAT, ROOT, COMM);
   /*
 * End of HPLAI_broadcast_AFLOAT
 */
}

#ifdef __cplusplus
}
#endif
