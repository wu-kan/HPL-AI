/*
 * Include files
 */
#include "hpl_ai.h"

template <typename T>
static void HPL_AI_atrsv_template(
	const enum HPL_ORDER ORDER,
	const enum HPL_UPLO UPLO,
	const enum HPL_TRANS TRANS,
	const enum HPL_DIAG DIAG,
	const int N,
	const T *A,
	const int LDA,
	T *X,
	const int INCX);

template <>
void HPL_AI_atrsv_template<double>(
	const enum HPL_ORDER ORDER,
	const enum HPL_UPLO UPLO,
	const enum HPL_TRANS TRANS,
	const enum HPL_DIAG DIAG,
	const int N,
	const double *A,
	const int LDA,
	double *X,
	const int INCX)
{
	HPL_dtrsv(ORDER, UPLO, TRANS, DIAG, N, A, LDA, X, INCX);
}

#if (HPL_CALL_CBLAS)
template <>
void HPL_AI_atrsv_template<float>(
	const enum HPL_ORDER ORDER,
	const enum HPL_UPLO UPLO,
	const enum HPL_TRANS TRANS,
	const enum HPL_DIAG DIAG,
	const int N,
	const float *A,
	const int LDA,
	float *X,
	const int INCX)
{
	cblas_strsv(ORDER, UPLO, TRANS, DIAG, N, A, LDA, X, INCX);
}
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
	void HPL_AI_atrsv(
		const enum HPL_ORDER ORDER,
		const enum HPL_UPLO UPLO,
		const enum HPL_TRANS TRANS,
		const enum HPL_DIAG DIAG,
		const int N,
		const HPL_AI_T_afloat *A,
		const int LDA,
		HPL_AI_T_afloat *X,
		const int INCX)
#else
void HPL_AI_atrsv(ORDER, UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
	const enum HPL_ORDER ORDER;
const enum HPL_UPLO UPLO;
const enum HPL_TRANS TRANS;
const enum HPL_DIAG DIAG;
const int N;
const HPL_AI_T_afloat *A;
const int LDA;
HPL_AI_T_afloat *X;
const int INCX;
#endif
	{
		HPL_AI_atrsv_template(ORDER, UPLO, TRANS, DIAG, N, A, LDA, X, INCX);
		/*
		 * End of HPL_AI_atrsv
 		 */
	}

#ifdef __cplusplus
}
#endif
