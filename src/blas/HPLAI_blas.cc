#include "hplai.hh"

template <>
int64_t blas::iamax<double>(
	int64_t n,
	double const *x,
	int64_t incx)
{
    //HPLAI_pabort( __LINE__, "blas::iamax<double>", "Use HPL_idamax" );
    //去掉上一行注释可以显示模板特化效果，下同
	return HPL_idamax(n, x, incx);
}

template <>
void blas::axpy<double, double>(
    int64_t n,
    blas::scalar_type<double, double> alpha,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy)
{
   HPL_daxpy(n, alpha, x, incx, y, incy);
}

template <>
void blas::copy<double, double>(
    int64_t n,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy)
{
    HPL_dcopy(n, x, incx, y, incy);
}

template <>
void blas::gemm<double, double, double>(
    blas::Layout layout,
    blas::Op transA,
    blas::Op transB,
    int64_t m,
    int64_t n,
    int64_t k,
    blas::scalar_type<double, double, double> alpha,
    double const *A,
    int64_t lda,
    double const *B,
    int64_t ldb,
    blas::scalar_type<double, double, double> beta,
    double *C,
    int64_t ldc)
{
    HPL_dgemm(
        layout == blas::Layout::RowMajor ? HplRowMajor : HplColumnMajor,
        transA == blas::Op::Trans ? HplTrans : HplNoTrans,
        transB == blas::Op::Trans ? HplTrans : HplNoTrans,
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

template <>
void blas::gemv<double, double, double>(
    blas::Layout layout,
    blas::Op trans,
    int64_t m,
    int64_t n,
    blas::scalar_type<double, double, double> alpha,
    double const *A,
    int64_t lda,
    double const *x,
    int64_t incx,
    blas::scalar_type<double, double, double> beta,
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

template <>
void blas::ger<double, double, double>(
    blas::Layout layout,
    int64_t m,
    int64_t n,
    blas::scalar_type<double, double, double> alpha,
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
        const_cast<double*>(y),
        incy,
        A,
        lda);
}

template <>
void blas::scal<double>(
    int64_t n,
    double alpha,
    double *x,
    int64_t incx)
{
    HPL_dscal(n, alpha, x, incx);
}

template <>
void blas::trsm<double, double>(
    blas::Layout layout,
    blas::Side side,
    blas::Uplo uplo,
    blas::Op trans,
    blas::Diag diag,
    int64_t m,
    int64_t n,
    blas::scalar_type<double, double> alpha,
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

template <>
void blas::trsv<double, double>(
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

#ifndef HPLAI_BLAS_GEN

template <>
int64_t blas::iamax<float>(
	int64_t n,
	float const *x,
	int64_t incx)
{
	return blas::iamax(n, x, incx);
}

template <>
void blas::axpy<float, float>(
    int64_t n,
    blas::scalar_type<float, float> alpha,
    float const *x,
    int64_t incx,
    float *y,
    int64_t incy)
{
   blas::axpy(n, alpha, x, incx, y, incy);
}

template <>
void blas::copy<float, float>(
    int64_t n,
    float const *x,
    int64_t incx,
    float *y,
    int64_t incy)
{
    blas::copy(n, x, incx, y, incy);
}

template <>
void blas::gemm<float, float, float>(
    blas::Layout layout,
    blas::Op transA,
    blas::Op transB,
    int64_t m,
    int64_t n,
    int64_t k,
    blas::scalar_type<float, float, float> alpha,
    float const *A,
    int64_t lda,
    float const *B,
    int64_t ldb,
    blas::scalar_type<float, float, float> beta,
    float *C,
    int64_t ldc)
{
    blas::gemm(
        layout,
        transA,
        transB,
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

template <>
void blas::gemv<float, float, float>(
    blas::Layout layout,
    blas::Op trans,
    int64_t m,
    int64_t n,
    blas::scalar_type<float, float, float> alpha,
    float const *A,
    int64_t lda,
    float const *x,
    int64_t incx,
    blas::scalar_type<float, float, float> beta,
    float *y,
    int64_t incy)
{
    blas::gemv(
        layout,
        trans,
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

template <>
void blas::ger<float, float, float>(
    blas::Layout layout,
    int64_t m,
    int64_t n,
    blas::scalar_type<float, float, float> alpha,
    float const *x,
    int64_t incx,
    float const *y,
    int64_t incy,
    float *A,
    int64_t lda)
{
    blas::ger(
        layout,
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

template <>
void blas::scal<float>(
    int64_t n,
    float alpha,
    float *x,
    int64_t incx)
{
    blas::scal(n, alpha, x, incx);
}

template <>
void blas::trsm<float, float>(
    blas::Layout layout,
    blas::Side side,
    blas::Uplo uplo,
    blas::Op trans,
    blas::Diag diag,
    int64_t m,
    int64_t n,
    blas::scalar_type<float, float> alpha,
    float const *A,
    int64_t lda,
    float *B,
    int64_t ldb)
{
    blas::trsm(
        layout,
        side,
        uplo,
        trans,
        diag,
        m,
        n,
        alpha,
        A,
        lda,
        B,
        ldb);
}

template <>
void blas::trsv<float, float>(
	blas::Layout layout,
	blas::Uplo uplo,
	blas::Op trans,
	blas::Diag diag,
	int64_t n,
	float const *A,
	int64_t lda,
	float *x,
	int64_t incx)
{
	blas::trsv(
		layout,
		uplo,
		trans,
		diag,
		n,
		A,
		lda,
		x,
		incx);
}

#else

template <>
void blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    blas::Op transA,
    blas::Op transB,
    int64_t m,
    int64_t n,
    int64_t k,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> alpha,
    HPLAI_T_AFLOAT const *A,
    int64_t lda,
    HPLAI_T_AFLOAT const *B,
    int64_t ldb,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> beta,
    HPLAI_T_AFLOAT *C,
    int64_t ldc)
{
    blas::gemm(
        layout,
        transA,
        transB,
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

template <>
void blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    blas::Side side,
    blas::Uplo uplo,
    blas::Op trans,
    blas::Diag diag,
    int64_t m,
    int64_t n,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> alpha,
    HPLAI_T_AFLOAT const *A,
    int64_t lda,
    HPLAI_T_AFLOAT *B,
    int64_t ldb)
{
    blas::trsm(
        layout,
        side,
        uplo,
        trans,
        diag,
        m,
        n,
        alpha,
        A,
        lda,
        B,
        ldb);
}

template <>
void blas::trsv<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
	blas::Layout layout,
	blas::Uplo uplo,
	blas::Op trans,
	blas::Diag diag,
	int64_t n,
	HPLAI_T_AFLOAT const *A,
	int64_t lda,
	HPLAI_T_AFLOAT *x,
	int64_t incx)
{
	blas::trsv(
		layout,
		uplo,
		trans,
		diag,
		n,
		A,
		lda,
		x,
		incx);
}

#endif

#ifdef __cplusplus
extern "C"
{
#endif

MPI_Datatype HPLAI_MPI_AFLOAT;

#ifdef STDC_HEADERS
    void HPLAI_blas_init(
        const int RANK,
        const int SIZE)
#else
void HPLAI_blas_init(RANK, SIZE)
    const int RANK,
    SIZE;
#endif
    {
        MPI_Type_contiguous(sizeof(HPLAI_T_AFLOAT), MPI_BYTE, &HPLAI_MPI_AFLOAT);
        MPI_Type_commit(&HPLAI_MPI_AFLOAT);
#ifdef HPL_CALL_VSIPL
        vsip_init((void *)0);
#endif
    }

    void HPLAI_blas_finalize()
    {
#ifdef HPL_CALL_VSIPL
        vsip_finalize((void *)0);
#endif
        MPI_Type_free(&HPLAI_MPI_AFLOAT);
    }

#ifdef __cplusplus
}
#endif