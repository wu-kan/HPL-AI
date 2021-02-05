/*
 * Include files
 */
#include "hpl_ai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
int HPL_AI_broadcast_afloat(
    void *BUFFER,
    const int COUNT,
    const int ROOT,
    MPI_Comm COMM)
#else
int HPL_AI_broadcast_afloat(BUFFER, COUNT, DTYPE, ROOT, COMM) void *BUFFER;
const int COUNT;
const int ROOT;
MPI_Comm COMM;
#endif
{
   return MPI_Bcast(BUFFER, 1LL * sizeof(HPL_AI_T_afloat) * COUNT, MPI_BYTE,
                    ROOT, COMM);
   /*
 * End of HPL_AI_broadcast_afloat
 */
}

#ifdef __cplusplus
}
#endif
