/*
 * Include files
 */
#include "hplai.h"

template <>
void blas::copy<double>(
    int64_t n,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy)
{
    HPL_dcopy(n, x, incx, y, incy);
}

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_acopy(
        const int N,
        const HPLAI_T_AFLOAT *X,
        const int INCX,
        HPLAI_T_AFLOAT *Y,
        const int INCY)
#else
void HPLAI_acopy(N, X, INCX, Y, INCY)
    const int N;
const HPLAI_T_AFLOAT *X;
const int INCX;
HPLAI_T_AFLOAT *Y;
const int INCY;
#endif
    {
        blas::copy(N, X, INCX, Y, INCY);
        /*
 * End of HPLAI_acopy
 */
    }

#ifdef __cplusplus
}
#endif
