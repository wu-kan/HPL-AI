/*
 * Include files
 */
#include "hplai.h"

template <typename T>
static void HPLAI_agemv_template(
    const enum HPL_ORDER ORDER,
    const enum HPL_TRANS TRANS,
    const int M,
    const int N,
    const T ALPHA,
    const T *A,
    const int LDA,
    const T *X,
    const int INCX,
    const T BETA,
    T *Y,
    const int INCY);

template <>
void HPLAI_agemv_template<double>(
    const enum HPL_ORDER ORDER,
    const enum HPL_TRANS TRANS,
    const int M,
    const int N,
    const double ALPHA,
    const double *A,
    const int LDA,
    const double *X,
    const int INCX,
    const double BETA,
    double *Y,
    const int INCY)
{
    HPL_dgemv(ORDER, TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
}

#ifdef HPL_CALL_CBLAS
template <>
void HPLAI_agemv_template<float>(
    const enum HPL_ORDER ORDER,
    const enum HPL_TRANS TRANS,
    const int M,
    const int N,
    const float ALPHA,
    const float *A,
    const int LDA,
    const float *X,
    const int INCX,
    const float BETA,
    float *Y,
    const int INCY)
{
    cblas_sgemv(ORDER, TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
}
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_agemv(
        const enum HPL_ORDER ORDER,
        const enum HPL_TRANS TRANS,
        const int M,
        const int N,
        const HPLAI_T_AFLOAT ALPHA,
        const HPLAI_T_AFLOAT *A,
        const int LDA,
        const HPLAI_T_AFLOAT *X,
        const int INCX,
        const HPLAI_T_AFLOAT BETA,
        HPLAI_T_AFLOAT *Y,
        const int INCY)
#else
void HPLAI_agemv(ORDER, TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
    const enum HPL_ORDER ORDER;
const enum HPL_TRANS TRANS;
const int M;
const int N;
const HPLAI_T_AFLOAT ALPHA;
const HPLAI_T_AFLOAT *A;
const int LDA;
const HPLAI_T_AFLOAT *X;
const int INCX;
const HPLAI_T_AFLOAT BETA;
HPLAI_T_AFLOAT *Y;
const int INCY;
#endif
    {
        HPLAI_agemv_template(ORDER, TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
        /*
 * End of HPLAI_agemv
 */
    }

#ifdef __cplusplus
}
#endif
