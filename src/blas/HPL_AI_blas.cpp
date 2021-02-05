#include "hpl_ai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPL_AI_blas_init(
        const int RANK,
        const int SIZE)
#else
void HPL_AI_blas_init(RANK, SIZE)
    const int RANK,
    SIZE;
#endif
    {
#ifdef HPL_CALL_VSIPL
        vsip_init((void *)0);
#endif
    }

    void HPL_AI_blas_finalize()
    {
#ifdef HPL_CALL_VSIPL
        vsip_finalize((void *)0);
#endif
    }

#ifdef __cplusplus
}
#endif