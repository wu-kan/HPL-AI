#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_blas_init(
        const int RANK,
        const int SIZE)
#else
void HPLAI_blas_init(RANK, SIZE)
    const int RANK,
    SIZE;
#endif
    {
#ifdef HPL_CALL_VSIPL
        vsip_init((void *)0);
#endif
    }

    void HPLAI_blas_finalize()
    {
#ifdef HPL_CALL_VSIPL
        vsip_finalize((void *)0);
#endif
    }

#ifdef __cplusplus
}
#endif