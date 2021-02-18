/*
 * Include files
 */
#include "hplai.h"

template <>
void blas::trsv<double>(
	blas::Layout layout,
	blas::Uplo uplo,
	blas::Op trans,
	blas::Diag diag,
	int64_t n,
	double const *A,
	int64_t lda,
	double *x,
	int64_t incx)
{
	HPL_dtrsv(
		layout == blas::Layout::RowMajor ? HplRowMajor : HplColumnMajor,
		uplo == blas::Uplo::Upper ? HplUpper : HplLower,
		trans == blas::Op::Trans ? HplTrans : HplNoTrans,
		diag == blas::Diag::Unit ? HplUnit : HplNonUnit,
		n,
		A,
		lda,
		x,
		incx);
}

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
		blas::trsv(
			ORDER == HPLAI_RowMajor ? blas::Layout::RowMajor : blas::Layout::ColMajor,
			UPLO == HPLAI_Upper ? blas::Uplo::Upper : blas::Uplo::Lower,
			TRANS == HPLAI_Trans ? blas::Op::Trans : blas::Op::NoTrans,
			DIAG == HPLAI_Unit ? blas::Diag::Unit : blas::Diag::NonUnit,
			N, A, LDA, X, INCX);
		/*
		 * End of HPLAI_atrsv
 		 */
	}

#ifdef __cplusplus
}
#endif
