/*
 * Include files
 */
#include "hplai.h"

template <typename T>
static void HPLAI_atrsv_template(
	const enum HPLAI_ORDER ORDER,
	const enum HPLAI_UPLO UPLO,
	const enum HPLAI_TRANS TRANS,
	const enum HPLAI_DIAG DIAG,
	const int N,
	const T *A,
	const int LDA,
	T *X,
	const int INCX);

template <>
void HPLAI_atrsv_template<double>(
	const enum HPLAI_ORDER ORDER,
	const enum HPLAI_UPLO UPLO,
	const enum HPLAI_TRANS TRANS,
	const enum HPLAI_DIAG DIAG,
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
void HPLAI_atrsv_template<float>(
	const enum HPLAI_ORDER ORDER,
	const enum HPLAI_UPLO UPLO,
	const enum HPLAI_TRANS TRANS,
	const enum HPLAI_DIAG DIAG,
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
	void HPLAI_atrsv(
		const enum HPLAI_ORDER ORDER,
		const enum HPLAI_UPLO UPLO,
		const enum HPLAI_TRANS TRANS,
		const enum HPLAI_DIAG DIAG,
		const int N,
		const HPLAI_T_AFLOAT *A,
		const int LDA,
		HPLAI_T_AFLOAT *X,
		const int INCX)
#else
void HPLAI_atrsv(ORDER, UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
	const enum HPLAI_ORDER ORDER;
const enum HPLAI_UPLO UPLO;
const enum HPLAI_TRANS TRANS;
const enum HPLAI_DIAG DIAG;
const int N;
const HPLAI_T_AFLOAT *A;
const int LDA;
HPLAI_T_AFLOAT *X;
const int INCX;
#endif
	{
		HPLAI_atrsv_template(ORDER, UPLO, TRANS, DIAG, N, A, LDA, X, INCX);
		/*
		 * End of HPLAI_atrsv
 		 */
	}

#ifdef __cplusplus
}
#endif
