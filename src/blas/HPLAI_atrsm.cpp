/*
 * Include files
 */
#include "hplai.h"

template <typename T>
static void HPLAI_atrsm_template(
    const enum HPLAI_ORDER ORDER,
    const enum HPLAI_SIDE SIDE,
    const enum HPLAI_UPLO UPLO,
    const enum HPLAI_TRANS TRANS,
    const enum HPLAI_DIAG DIAG,
    const int M,
    const int N,
    const T ALPHA,
    const T *A,
    const int LDA,
    T *B,
    const int LDB);

template <>
void HPLAI_atrsm_template<double>(
    const enum HPLAI_ORDER ORDER,
    const enum HPLAI_SIDE SIDE,
    const enum HPLAI_UPLO UPLO,
    const enum HPLAI_TRANS TRANS,
    const enum HPLAI_DIAG DIAG,
    const int M,
    const int N,
    const double ALPHA,
    const double *A,
    const int LDA,
    double *B,
    const int LDB)
{
    HPL_dtrsm(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB);
}

#ifdef HPL_CALL_CBLAS
template <>
void HPLAI_atrsm_template<float>(
    const enum HPLAI_ORDER ORDER,
    const enum HPLAI_SIDE SIDE,
    const enum HPLAI_UPLO UPLO,
    const enum HPLAI_TRANS TRANS,
    const enum HPLAI_DIAG DIAG,
    const int M,
    const int N,
    const float ALPHA,
    const float *A,
    const int LDA,
    float *B,
    const int LDB)
{
    cblas_strsm(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB);
}
#endif

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
        HPLAI_atrsm_template(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB);
        /*
 * End of HPLAI_atrsm
 */
    }

#ifdef __cplusplus
}
#endif
