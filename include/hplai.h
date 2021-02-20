#ifndef HPLAI_H
#define HPLAI_H

#ifdef __cplusplus
extern "C"
{
#endif
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl.h"

#include "hplai_misc.h"

#include "hplai_auxil.h"
#include "hplai_gesv.h"

#include "hplai_pauxil.h"
#include "hplai_panel.h"
#include "hplai_pfact.h"

#include "hplai_pmatgen.h"
#include "hplai_ptest.h"

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

//use blaspp https://bitbucket.org/icl/blaspp/src/master/
#include <blas.hh>
#include "hplai_pgesv.hh"

template <>
int64_t blas::iamax<double>(
	int64_t n,
	double const *x,
	int64_t incx);

template <>
void blas::axpy<double, double>(
    int64_t n,
    double alpha,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy);

template <>
void blas::copy<double>(
    int64_t n,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy);

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
    int64_t ldc);

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
    int64_t incy);

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
    int64_t lda);

template <>
void blas::scal<double>(
    int64_t n,
    double alpha,
    double *x,
    int64_t incx);
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
    int64_t ldb);

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
	int64_t incx);

#endif

#endif
/*
 * End of hplai.h
 */
