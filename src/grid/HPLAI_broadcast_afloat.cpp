/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
int HPLAI_broadcast_afloat(
    void *BUFFER,
    const int COUNT,
    const int ROOT,
    MPI_Comm COMM)
#else
int HPLAI_broadcast_afloat(BUFFER, COUNT, DTYPE, ROOT, COMM) void *BUFFER;
const int COUNT;
const int ROOT;
MPI_Comm COMM;
#endif
{
   return MPI_Bcast(BUFFER, 1LL * sizeof(HPLAI_T_AFLOAT) * COUNT, MPI_BYTE,
                    ROOT, COMM);
   /*
 * End of HPLAI_broadcast_afloat
 */
}

#ifdef __cplusplus
}
#endif
