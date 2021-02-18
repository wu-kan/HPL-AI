/*
 * Include files
 */
#include "hplai.h"

template <>
void blas::trsm<double, double>(
    blas::Layout layout,
    blas::Side side,
    blas::Uplo uplo,
    blas::Op trans,
    blas::Diag diag,
    int64_t m,
    int64_t n,
    double alpha,
    double const *A,
    int64_t lda,
    double *B,
    int64_t ldb)
{
    HPL_dtrsm(
        layout == blas::Layout::RowMajor ? HplRowMajor : HplColumnMajor,
        side == blas::Side::Left ? HplLeft : HplRight,
        uplo == blas::Uplo::Upper ? HplUpper : HplLower,
        trans == blas::Op::Trans ? HplTrans : HplNoTrans,
        diag == blas::Diag::Unit ? HplUnit : HplNonUnit,
        m,
        n,
        alpha,
        A,
        lda,
        B,
        ldb);
}

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_atrsm(
        const enum HPLAI_ORDER ORDER,
        const enum HPLAI_SIDE SIDE,
        const enum HPLAI_UPLO UPLO,
        const enum HPLAI_TRANS TRANS,
        const enum HPLAI_DIAG DIAG,
        const int M,
        const int N,
        const HPLAI_T_AFLOAT ALPHA,
        const HPLAI_T_AFLOAT *A,
        const int LDA,
        HPLAI_T_AFLOAT *B,
        const int LDB)
#else
void HPL_dtrsm(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB)
    const enum HPLAI_ORDER ORDER;
const enum HPLAI_SIDE SIDE;
const enum HPLAI_UPLO UPLO;
const enum HPLAI_TRANS TRANS;
const enum HPLAI_DIAG DIAG;
const int M;
const int N;
const HPLAI_T_AFLOAT ALPHA;
const HPLAI_T_AFLOAT *A;
const int LDA;
HPLAI_T_AFLOAT *B;
const int LDB;
#endif
    {
        blas::trsm(
            ORDER == HPLAI_RowMajor ? blas::Layout::RowMajor : blas::Layout::ColMajor,
            SIDE == HPLAI_Left ? blas::Side::Left : blas::Side::Right,
            UPLO == HPLAI_Upper ? blas::Uplo::Upper : blas::Uplo::Lower,
            TRANS == HPLAI_Trans ? blas::Op::Trans : blas::Op::NoTrans,
            DIAG == HPLAI_Unit ? blas::Diag::Unit : blas::Diag::NonUnit,
            M, N, ALPHA, A, LDA, B, LDB);
        /*
 * End of HPLAI_atrsm
 */
    }

#ifdef __cplusplus
}
#endif
