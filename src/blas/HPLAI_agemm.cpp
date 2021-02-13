/*
 * Include files
 */
#include "hplai.h"

template <typename T>
static void HPLAI_agemm_template(
    const enum HPL_ORDER ORDER,
    const enum HPL_TRANS TRANSA,
    const enum HPL_TRANS TRANSB,
    const int M,
    const int N,
    const int K,
    const T ALPHA,
    const T *A,
    const int LDA,
    const T *B,
    const int LDB,
    const T BETA,
    T *C,
    const int LDC);

template <>
void HPLAI_agemm_template<double>(
    const enum HPL_ORDER ORDER,
    const enum HPL_TRANS TRANSA,
    const enum HPL_TRANS TRANSB,
    const int M,
    const int N,
    const int K,
    const double ALPHA,
    const double *A,
    const int LDA,
    const double *B,
    const int LDB,
    const double BETA,
    double *C,
    const int LDC)
{
    HPL_dgemm(ORDER, TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
              BETA, C, LDC);
}

#ifdef HPL_CALL_CBLAS
template <>
void HPLAI_agemm_template<float>(
    const enum HPL_ORDER ORDER,
    const enum HPL_TRANS TRANSA,
    const enum HPL_TRANS TRANSB,
    const int M,
    const int N,
    const int K,
    const float ALPHA,
    const float *A,
    const int LDA,
    const float *B,
    const int LDB,
    const float BETA,
    float *C,
    const int LDC)
{
    cblas_sgemm(ORDER, TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
                BETA, C, LDC);
}
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_agemm(
        const enum HPL_ORDER ORDER,
        const enum HPL_TRANS TRANSA,
        const enum HPL_TRANS TRANSB,
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
    const enum HPL_ORDER ORDER;
const enum HPL_TRANS TRANSA;
const enum HPL_TRANS TRANSB;
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
        HPLAI_agemm_template(ORDER, TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
        /*
 * End of HPLAI_agemm
 */
    }

#ifdef __cplusplus
}
#endif
