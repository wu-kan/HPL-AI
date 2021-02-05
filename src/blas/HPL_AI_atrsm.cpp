/*
 * Include files
 */
#include "hpl_ai.h"

template <typename T>
static void HPL_AI_atrsm_template(
    const enum HPL_ORDER ORDER,
    const enum HPL_SIDE SIDE,
    const enum HPL_UPLO UPLO,
    const enum HPL_TRANS TRANS,
    const enum HPL_DIAG DIAG,
    const int M,
    const int N,
    const T ALPHA,
    const T *A,
    const int LDA,
    T *B,
    const int LDB);

template <>
void HPL_AI_atrsm_template<double>(
    const enum HPL_ORDER ORDER,
    const enum HPL_SIDE SIDE,
    const enum HPL_UPLO UPLO,
    const enum HPL_TRANS TRANS,
    const enum HPL_DIAG DIAG,
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
void HPL_AI_atrsm_template<float>(
    const enum HPL_ORDER ORDER,
    const enum HPL_SIDE SIDE,
    const enum HPL_UPLO UPLO,
    const enum HPL_TRANS TRANS,
    const enum HPL_DIAG DIAG,
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
    void HPL_AI_atrsm(
        const enum HPL_ORDER ORDER,
        const enum HPL_SIDE SIDE,
        const enum HPL_UPLO UPLO,
        const enum HPL_TRANS TRANS,
        const enum HPL_DIAG DIAG,
        const int M,
        const int N,
        const HPL_AI_T_afloat ALPHA,
        const HPL_AI_T_afloat *A,
        const int LDA,
        HPL_AI_T_afloat *B,
        const int LDB)
#else
void HPL_dtrsm(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB)
    const enum HPL_ORDER ORDER;
const enum HPL_SIDE SIDE;
const enum HPL_UPLO UPLO;
const enum HPL_TRANS TRANS;
const enum HPL_DIAG DIAG;
const int M;
const int N;
const HPL_AI_T_afloat ALPHA;
const HPL_AI_T_afloat *A;
const int LDA;
HPL_AI_T_afloat *B;
const int LDB;
#endif
    {
        HPL_AI_atrsm_template(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB);
        /*
 * End of HPL_AI_atrsm
 */
    }

#ifdef __cplusplus
}
#endif
