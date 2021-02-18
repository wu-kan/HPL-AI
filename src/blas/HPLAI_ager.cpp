/*
 * Include files
 */
#include "hplai.h"

template <>
void blas::ger<double, double, double>(
    blas::Layout layout,
    int64_t m,
    int64_t n,
    double alpha,
    double const *x,
    int64_t incx,
    double const *y,
    int64_t incy,
    double *A,
    int64_t lda)
{
    HPL_dger(
        layout == blas::Layout::RowMajor ? HplRowMajor : HplColumnMajor,
        m,
        n,
        alpha,
        x,
        incx,
        y,
        incy,
        A,
        lda);
}

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_ager(
        const enum HPLAI_ORDER ORDER,
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
    const enum HPLAI_ORDER ORDER;
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
        blas::ger(
            ORDER == HPLAI_RowMajor ? blas::Layout::RowMajor : blas::Layout::ColMajor,
            M, N, ALPHA, X, INCX, Y, INCY, A, LDA);
        /*
 * End of HPLAI_ager
 */
    }

#ifdef __cplusplus
}
#endif
