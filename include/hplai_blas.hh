/*
 * MIT License
 * 
 * Copyright (c) 2021 WuK
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef HPLAI_BLAS_HH
#define HPLAI_BLAS_HH

#include "hplai_misc.hh"

#ifdef __cplusplus

template <>
int64_t blas::iamax<double>(
    int64_t n,
    double const *x,
    int64_t incx);

template <>
int64_t blas::iamax<HPLAI_T_AFLOAT>(
    int64_t n,
    HPLAI_T_AFLOAT const *x,
    int64_t incx);

template <>
void blas::axpy<double, double>(
    int64_t n,
    blas::scalar_type<double, double> alpha,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy);

template <>
void blas::axpy<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    int64_t n,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> alpha,
    HPLAI_T_AFLOAT const *x,
    int64_t incx,
    HPLAI_T_AFLOAT *y,
    int64_t incy);

template <>
void blas::copy<double, double>(
    int64_t n,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy);

template <>
void blas::copy<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    int64_t n,
    HPLAI_T_AFLOAT const *x,
    int64_t incx,
    HPLAI_T_AFLOAT *y,
    int64_t incy);

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
    int64_t ldc);

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
    int64_t ldc);

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
    int64_t incy);

template <>
void blas::gemv<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    blas::Op trans,
    int64_t m,
    int64_t n,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> alpha,
    HPLAI_T_AFLOAT const *A,
    int64_t lda,
    HPLAI_T_AFLOAT const *x,
    int64_t incx,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> beta,
    HPLAI_T_AFLOAT *y,
    int64_t incy);

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
    int64_t lda);

template <>
void blas::ger<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    int64_t m,
    int64_t n,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> alpha,
    HPLAI_T_AFLOAT const *x,
    int64_t incx,
    HPLAI_T_AFLOAT const *y,
    int64_t incy,
    HPLAI_T_AFLOAT *A,
    int64_t lda);

template <>
void blas::scal<double>(
    int64_t n,
    double alpha,
    double *x,
    int64_t incx);

template <>
void blas::scal<HPLAI_T_AFLOAT>(
    int64_t n,
    HPLAI_T_AFLOAT alpha,
    HPLAI_T_AFLOAT *x,
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
    blas::scalar_type<double, double> alpha,
    double const *A,
    int64_t lda,
    double *B,
    int64_t ldb);

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
    int64_t ldb);

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
    int64_t incx);

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
    int64_t incx);

#endif

#ifdef __cplusplus
extern "C"
{
#endif

    extern MPI_Datatype HPLAI_MPI_AFLOAT;

    void HPLAI_blas_init
        STDC_ARGS((
            const int,
            const int));

    void HPLAI_blas_finalize();

#ifdef __cplusplus
}
#endif

#endif
/*
 * End of hplai_blas.hh
 */
