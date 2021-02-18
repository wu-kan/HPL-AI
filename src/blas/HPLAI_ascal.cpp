/*
 * Include files
 */
#include "hplai.h"

template <>
void blas::scal<double>(
    int64_t n,
    double alpha,
    double *x,
    int64_t incx)
{
    HPL_dscal(n, alpha, x, incx);
}

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_ascal(
        const int N,
        const HPLAI_T_AFLOAT ALPHA,
        HPLAI_T_AFLOAT *X,
        const int INCX)
#else
void HPLAI_ascal(N, ALPHA, X, INCX)
    const int N;
const HPLAI_T_AFLOAT ALPHA;
HPLAI_T_AFLOAT *X;
const int INCX;
#endif
    {
        blas::scal(N, ALPHA, X, INCX);
        /*
 * End of HPLAI_ascal
 */
    }

#ifdef __cplusplus
}
#endif
