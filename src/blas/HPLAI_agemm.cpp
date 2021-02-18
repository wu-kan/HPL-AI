/*
 * Include files
 */
#include "hplai.h"

template <>
void blas::gemm<double, double, double>(
    blas::Layout layout,
    blas::Op transA,
    blas::Op transB,
    int64_t m,
    int64_t n,
    int64_t k,
    double alpha,
    double const *A,
    int64_t lda,
    double const *B,
    int64_t ldb,
    double beta,
    double *C,
    int64_t ldc)
{
    HPL_dgemm(
        layout == blas::Layout::RowMajor ? HPLAI_RowMajor : HPLAI_ColumnMajor,
        transA == blas::Op::Trans ? HPLAI_Trans : HPLAI_NoTrans,
        transB == blas::Op::Trans ? HPLAI_Trans : HPLAI_NoTrans,
        m,
        n,
        k,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc);
}

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_agemm(
        const enum HPLAI_ORDER ORDER,
        const enum HPLAI_TRANS TRANSA,
        const enum HPLAI_TRANS TRANSB,
        const int M,
        const int N,
        const int K,
        const HPLAI_T_AFLOAT ALPHA,
        const HPLAI_T_AFLOAT *A,
        const int LDA,
        const HPLAI_T_AFLOAT *B,
        const int LDB,
        const HPLAI_T_AFLOAT BETA,
        HPLAI_T_AFLOAT *C,
        const int LDC)
#else
void HPLAI_agemm(ORDER, TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    const enum HPLAI_ORDER ORDER;
const enum HPLAI_TRANS TRANSA;
const enum HPLAI_TRANS TRANSB;
const int M;
const int N;
const int K;
const HPLAI_T_AFLOAT ALPHA;
const HPLAI_T_AFLOAT *A;
const int LDA;
const HPLAI_T_AFLOAT *B;
const int LDB;
const HPLAI_T_AFLOAT BETA;
HPLAI_T_AFLOAT *C;
const int LDC;
#endif
    {
        blas::gemm(
            ORDER == HPLAI_RowMajor ? blas::Layout::RowMajor : blas::Layout::ColMajor,
            TRANSA == HPLAI_Trans ? blas::Op::Trans : blas::Op::NoTrans,
            TRANSB == HPLAI_Trans ? blas::Op::Trans : blas::Op::NoTrans,
            M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
        /*
 * End of HPLAI_agemm
 */
    }

#ifdef __cplusplus
}
#endif
