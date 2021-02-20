#include "hplai.h"

template <>
void blas::axpy<double, double>(
    int64_t n,
    double alpha,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy)
{
   HPL_daxpy(n, alpha, x, incx, y, incy);
}

template <>
void blas::copy<double>(
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
    double alpha,
    double const *A,
    int64_t lda,
    double const *x,
    int64_t incx,
    double beta,
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

MPI_Datatype HPLAI_MPI_AFLOAT;

#ifdef STDC_HEADERS
    void HPLAI_init(
        const int RANK,
        const int SIZE)
#else
void HPLAI_init(RANK, SIZE)
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

    void HPLAI_finalize()
    {
#ifdef HPL_CALL_VSIPL
        vsip_finalize((void *)0);
#endif
        MPI_Type_free(&HPLAI_MPI_AFLOAT);
    }

#ifdef __cplusplus
}
#endif