/*
 * Include files
 */
#include "hplai.h"

template <typename T>
static void HPLAI_ager_template(
    const enum HPL_ORDER ORDER,
    const int M,
    const int N,
    const T ALPHA,
    const T *X,
    const int INCX,
    T *Y,
    const int INCY,
    T *A,
    const int LDA);

template <>
void HPLAI_ager_template<double>(
    const enum HPL_ORDER ORDER,
    const int M,
    const int N,
    const double ALPHA,
    const double *X,
    const int INCX,
    double *Y,
    const int INCY,
    double *A,
    const int LDA)
{
    HPL_dger(ORDER, M, N, ALPHA, X, INCX, Y, INCY, A, LDA);
}

#ifdef HPL_CALL_CBLAS
template <>
void HPLAI_ager_template<float>(
    const enum HPL_ORDER ORDER,
    const int M,
    const int N,
    const float ALPHA,
    const float *X,
    const int INCX,
    float *Y,
    const int INCY,
    float *A,
    const int LDA)
{
    cblas_sger(ORDER, M, N, ALPHA, X, INCX, Y, INCY, A, LDA);
}
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_ager(
        const enum HPL_ORDER ORDER,
        const int M,
        const int N,
        const HPLAI_T_AFLOAT ALPHA,
        const HPLAI_T_AFLOAT *X,
        const int INCX,
        HPLAI_T_AFLOAT *Y,
        const int INCY,
        HPLAI_T_AFLOAT *A,
        const int LDA)
#else
void HPLAI_ager(ORDER, M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
    const enum HPL_ORDER ORDER;
const int M;
const int N;
const HPLAI_T_AFLOAT ALPHA;
const HPLAI_T_AFLOAT *X;
const int INCX;
HPLAI_T_AFLOAT *Y;
const int INCY;
HPLAI_T_AFLOAT *A;
const int LDA;
#endif
    {
        HPLAI_ager_template(ORDER, M, N, ALPHA, X, INCX, Y, INCY, A, LDA);
        /*
 * End of HPLAI_ager
 */
    }

#ifdef __cplusplus
}
#endif
