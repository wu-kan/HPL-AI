/*
 * Include files
 */
#include "hplai.h"

template <>
void blas::gemv<double, double, double>(
    blas::Layout layout,
    blas::Op trans,
    int64_t m,
    int64_t n,
    double alpha,
    double const *A,
    int64_t lda,
    double const *x,
    int64_t incx,
    double beta,
    double *y,
    int64_t incy)
{
    HPL_dgemv(
        layout == blas::Layout::RowMajor ? HplRowMajor : HplColumnMajor,
        trans == blas::Op::Trans ? HplTrans : HplNoTrans,
        m,
        n,
        alpha,
        A,
        lda,
        x,
        incx,
        beta,
        y,
        incy);
}

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_agemv(
        const enum HPLAI_ORDER ORDER,
        const enum HPLAI_TRANS TRANS,
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
    const enum HPLAI_ORDER ORDER;
const enum HPLAI_TRANS TRANS;
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
        blas::gemv(
            ORDER == HPLAI_RowMajor ? blas::Layout::RowMajor : blas::Layout::ColMajor,
            TRANS == HPLAI_Trans ? blas::Op::Trans : blas::Op::NoTrans,
            M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY);
        /*
 * End of HPLAI_agemv
 */
    }

#ifdef __cplusplus
}
#endif
