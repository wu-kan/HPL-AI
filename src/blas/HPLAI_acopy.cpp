/*
 * Include files
 */
#include "hplai.h"

template <typename T>
static void HPLAI_acopy_template(
    const int N,
    const T *X,
    const int INCX,
    T *Y,
    const int INCY);

template <>
void HPLAI_acopy_template<double>(
    const int N,
    const double *X,
    const int INCX,
    double *Y,
    const int INCY)
{
    HPL_dcopy(N, X, INCX, Y, INCY);
}

#ifdef HPL_CALL_CBLAS
template <>
void HPLAI_acopy_template<float>(
    const int N,
    const float *X,
    const int INCX,
    float *Y,
    const int INCY)
{
    cblas_scopy(N, X, INCX, Y, INCY);
}
#endif

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
        HPLAI_acopy_template(N, X, INCX, Y, INCY);
        /*
 * End of HPLAI_acopy
 */
    }

#ifdef __cplusplus
}
#endif
