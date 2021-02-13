/*
 * Include files
 */
#include "hplai.h"

template <typename T>
static void HPLAI_ascal_template(
    const int N,
    const T ALPHA,
    T *X,
    const int INCX);

template <>
void HPLAI_ascal_template<double>(
    const int N,
    const double ALPHA,
    double *X,
    const int INCX)
{
    HPL_dscal(N, ALPHA, X, INCX);
}

#ifdef HPL_CALL_CBLAS
template <>
void HPLAI_ascal_template<float>(
    const int N,
    const float ALPHA,
    float *X,
    const int INCX)
{
    cblas_sscal(N, ALPHA, X, INCX);
}
#endif

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
        HPLAI_ascal_template(N, ALPHA, X, INCX);
        /*
 * End of HPLAI_ascal
 */
    }

#ifdef __cplusplus
}
#endif
