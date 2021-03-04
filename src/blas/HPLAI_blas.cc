#include "hplai.hh"

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

template <>
int64_t blas::iamax<double>(
	int64_t n,
	double const *x,
	int64_t incx)
{
	return HPL_idamax(n, x, incx);
}

template <>
int64_t blas::iamax<float>(
	int64_t n,
	float const *x,
	int64_t incx)
{
	return blas::iamax(n, x, incx);
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

#ifndef HPLAI_GEN_BLASPP_GEMM

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

#else

template<typename T>
static  void HPLAI_gemmNN
(
   const int64_t                  M,
   const int64_t                  N,
   const int64_t                  K,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   const T               * B,
   const int64_t                  LDB,
   const T               BETA,
   T                     * C,
   const int64_t                  LDC
)
{
   register T            t0;
   int64_t                        i, iail, iblj, icij, j, jal, jbj, jcj, l;

   for( j = 0, jbj = 0, jcj  = 0; j < N; j++, jbj += LDB, jcj += LDC )
   {
      blas::scal<T>( M, BETA, C+jcj, 1 );
      for( l = 0, jal = 0, iblj = jbj; l < K; l++, jal += LDA, iblj += 1 )
      {
         t0 = ALPHA * B[iblj];
         for( i = 0, iail = jal, icij = jcj; i < M; i++, iail += 1, icij += 1 )
         { C[icij] += A[iail] * t0; }
      }
   }
}

template<typename T>
static  void HPLAI_gemmNT
(
   const int64_t                  M,
   const int64_t                  N,
   const int64_t                  K,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   const T               * B,
   const int64_t                  LDB,
   const T               BETA,
   T                     * C,
   const int64_t                  LDC
)
{
   register T            t0;
   int64_t                        i, iail, ibj, ibjl, icij, j, jal, jcj, l;

   for( j = 0, ibj  = 0, jcj  = 0; j < N; j++, ibj += 1, jcj += LDC )
   {
      blas::scal<T>( M, BETA, C+jcj, 1 );
      for( l = 0, jal = 0, ibjl = ibj; l < K; l++, jal += LDA, ibjl += LDB )
      {
         t0 = ALPHA * B[ibjl];
         for( i = 0, iail = jal, icij = jcj; i < M; i++, iail += 1, icij += 1 )
         { C[icij] += A[iail] * t0; }
      }
   }
}

template<typename T>
static  void HPLAI_gemmTN
(
   const int64_t                  M,
   const int64_t                  N,
   const int64_t                  K,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   const T               * B,
   const int64_t                  LDB,
   const T               BETA,
   T                     * C,
   const int64_t                  LDC
)
{
   register T            t0;
   int64_t                        i, iai, iail, iblj, icij, j, jbj, jcj, l;

   for( j = 0, jbj = 0, jcj = 0; j < N; j++, jbj += LDB, jcj += LDC )
   {
      for( i = 0, icij = jcj, iai = 0; i < M; i++, icij += 1, iai += LDA )
      {
         t0 = HPLAI_rzero;
         for( l = 0, iail = iai, iblj = jbj; l < K; l++, iail += 1, iblj += 1 )
         { t0 += A[iail] * B[iblj]; }
         if( BETA == HPLAI_rzero ) C[icij]  = HPLAI_rzero;
         else                    C[icij] *= BETA;
         C[icij] += ALPHA * t0;
      }
   }
}

template<typename T>
static  void HPLAI_gemmTT
(
   const int64_t                  M,
   const int64_t                  N,
   const int64_t                  K,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   const T               * B,
   const int64_t                  LDB,
   const T               BETA,
   T                     * C,
   const int64_t                  LDC
)
{
   register T            t0;
   int64_t                        i, iali, ibj, ibjl, icij, j, jai, jcj, l;

   for( j = 0, ibj = 0, jcj  = 0; j < N; j++, ibj += 1, jcj += LDC )
   {
      for( i = 0, icij = jcj, jai = 0; i < M; i++, icij += 1, jai += LDA )
      {
         t0 = HPLAI_rzero;
         for( l = 0,      iali  = jai, ibjl  = ibj;
              l < K; l++, iali += 1,   ibjl += LDB ) t0 += A[iali] * B[ibjl];
         if( BETA == HPLAI_rzero ) C[icij]  = HPLAI_rzero;
         else                    C[icij] *= BETA;
         C[icij] += ALPHA * t0;
      }
   }
}

template <>
void blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    blas::Op TRANSA,
    blas::Op TRANSB,
    int64_t M,
    int64_t N,
    int64_t K,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> ALPHA,
    HPLAI_T_AFLOAT const *A,
    int64_t LDA,
    HPLAI_T_AFLOAT const *B,
    int64_t LDB,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> BETA,
    HPLAI_T_AFLOAT *C,
    int64_t LDC)
{
   //HPLAI_pabort( __LINE__, "blas::gemm", "Use HPLAI_GEN_BLASPP_GEMM" );
   if(layout != blas::Layout::ColMajor)
   {
      blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
         blas::Layout::ColMajor,
         TRANSB, TRANSA, N, M, K, ALPHA, B, LDB, A, LDA, BETA, C, LDC );
      return;
   }
   int64_t                        i, j;

   if( ( M == 0 ) || ( N == 0 ) ||
       ( ( ( ALPHA == HPLAI_rzero ) || ( K == 0 ) ) &&
         ( BETA == HPLAI_rone ) ) ) return;

   if( ALPHA == HPLAI_rzero )
   {
      for( j = 0; j < N; j++ )
      {  for( i = 0; i < M; i++ ) *(C+i+j*LDC) = HPLAI_rzero; }
      return;
   }

   if( TRANSB == blas::Op::NoTrans )
   {
      if( TRANSA == blas::Op::NoTrans )
      { HPLAI_gemmNN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
      else
      { HPLAI_gemmTN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
   }
   else
   {
      if( TRANSA == blas::Op::NoTrans )
      { HPLAI_gemmNT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
      else
      { HPLAI_gemmTT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
   }
}

#endif

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
void blas::scal<double>(
    int64_t n,
    double alpha,
    double *x,
    int64_t incx)
{
    HPL_dscal(n, alpha, x, incx);
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

#ifndef HPLAI_GEN_BLASPP_TRSM

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

#else

template<typename T>
static  void HPLAI_trsmLLNN
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   int64_t                        i, iaik, ibij, ibkj, j, jak, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, ibij= jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = 0, jak  = 0, ibkj = jbj; k < M; k++, jak += LDA, ibkj += 1 )
      {
         B[ibkj] /= A[k+jak];
         for( i = k+1,    iaik  = k+1+jak, ibij  = k+1+jbj;
              i < M; i++, iaik +=1,        ibij += 1 )
         { B[ibij] -= B[ibkj] * A[iaik]; }
      }
   }
}

template<typename T>
static  void HPLAI_trsmLLNU
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   int64_t                        i, iaik, ibij, ibkj, j, jak, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, ibij= jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = 0, jak  = 0, ibkj = jbj; k < M; k++, jak += LDA, ibkj += 1 )
      {
         for( i = k+1,    iaik  = k+1+jak, ibij  = k+1+jbj;
              i < M; i++, iaik +=1,        ibij += 1 )
         { B[ibij] -= B[ibkj] * A[iaik]; }
      }
   }
}

template<typename T>
static  void HPLAI_trsmLLTN
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   register T            t0;
   int64_t                        i, iaki, ibij, ibkj, j, jai, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = M-1,     jai  = (M-1)*LDA, ibij  = M-1+jbj;
           i >= 0; i--, jai -= LDA,       ibij -= 1 )
      {
         t0 = ALPHA * B[ibij];
         for( k = i+1,    iaki  = i+1+jai, ibkj  = i+1+jbj;
              k < M; k++, iaki += 1,       ibkj += 1 )
         { t0 -= A[iaki] * B[ibkj]; }
         t0 /= A[i+jai];
         B[ibij] = t0;
      }
   }
}

template<typename T>
static  void HPLAI_trsmLLTU
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   register T            t0;
   int64_t                        i, iaki, ibij, ibkj, j, jai, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = M-1,     jai  = (M-1)*LDA, ibij  = M-1+jbj;
           i >= 0; i--, jai -= LDA,       ibij -= 1 )
      {
         t0 = ALPHA * B[ibij];
         for( k = i+1,    iaki  = i+1+jai, ibkj  = i+1+jbj;
              k < M; k++, iaki += 1,       ibkj += 1 )
         { t0 -= A[iaki] * B[ibkj]; }
         B[ibij] = t0;
      }
   }
}

template<typename T>
static  void HPLAI_trsmLUNN
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   int64_t                        i, iaik, ibij, ibkj, j, jak, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = M-1,     jak  = (M-1)*LDA, ibkj  = M-1+jbj;
           k >= 0; k--, jak -= LDA,       ibkj -= 1 )
      {
         B[ibkj] /= A[k+jak];
         for( i = 0,      iaik  = jak, ibij  = jbj;
              i < k; i++, iaik += 1,   ibij += 1 )
         { B[ibij] -= B[ibkj] * A[iaik]; }
      }
   }
}

template<typename T>
static  void HPLAI_trsmLUNU
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   int64_t                        i, iaik, ibij, ibkj, j, jak, jbj, k;

   for( j = 0, jbj = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = M-1,     jak  = (M-1)*LDA, ibkj  = M-1+jbj;
           k >= 0; k--, jak -= LDA,       ibkj -= 1 )
      {
         for( i = 0,      iaik  = jak, ibij  = jbj;
              i < k; i++, iaik += 1,   ibij += 1 )
         { B[ibij] -= B[ibkj] * A[iaik]; }
      }
   }
}

template<typename T>
static  void HPLAI_trsmLUTN
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   int64_t                        i, iaki, ibij, ibkj, j, jai, jbj, k;
   register T            t0;

   for( j = 0, jbj  = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, jai  = 0, ibij = jbj; i < M; i++, jai += LDA, ibij += 1 )
      {
         t0 = ALPHA * B[ibij];
         for( k = 0, iaki = jai, ibkj = jbj; k < i; k++, iaki += 1, ibkj += 1 )
         { t0 -= A[iaki] * B[ibkj]; }
         t0 /= A[i+jai];
         B[ibij] = t0;
      }
   }
}

template<typename T>
static  void HPLAI_trsmLUTU
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   register T            t0;
   int64_t                        i, iaki, ibij, ibkj, j, jai, jbj, k;

   for( j = 0, jbj  = 0; j < N; j++, jbj += LDB )
   {
      for( i = 0, jai  = 0, ibij = jbj; i < M; i++, jai += LDA, ibij += 1 )
      {
         t0 = ALPHA * B[ibij];
         for( k = 0, iaki = jai, ibkj = jbj; k < i; k++, iaki += 1, ibkj += 1 )
         { t0 -= A[iaki] * B[ibkj]; }
         B[ibij] = t0;
      }
   }
}

template<typename T>
static  void HPLAI_trsmRLNN
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   int64_t                        i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

   for( j = N-1,      jaj  = (N-1)*LDA, jbj  = (N-1)*LDB;
        j >= 0;  j--, jaj -= LDA,       jbj -= LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = j+1,    iakj  = j+1+jaj, jbk  = (j+1)*LDB;
           k < N; k++, iakj += 1,       jbk += LDB )
      {
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= A[iakj] * B[ibik]; }
      }
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] /= A[j+jaj]; }
   }
}

template<typename T>
static  void HPLAI_trsmRLNU
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   int64_t                        i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

   for( j = N-1,      jaj  = (N-1)*LDA, jbj  = (N-1)*LDB;
        j >= 0;  j--, jaj -= LDA,       jbj -= LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = j+1,    iakj  = j+1+jaj, jbk  = (j+1)*LDB;
           k < N; k++, iakj += 1,       jbk += LDB )
      {
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= A[iakj] * B[ibik]; }
      }
   }
}

template<typename T>
static  void HPLAI_trsmRLTN
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   register T            t0;
   int64_t                        i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

   for( k = 0, jak = 0, jbk = 0; k < N; k++, jak += LDA, jbk += LDB )
   {
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] /= A[k+jak]; }
      for( j = k+1,    iajk  = (k+1)+jak, jbj  = (k+1)*LDB;
           j < N; j++, iajk += 1,         jbj += LDB )
      {
         t0 = A[iajk];
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= t0 * B[ibik]; }
      }
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] *= ALPHA; }
   }
}

template<typename T>
static  void HPLAI_trsmRLTU
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   register T            t0;
   int64_t                        i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

   for( k = 0, jak = 0, jbk = 0; k < N; k++, jak += LDA, jbk += LDB )
   {
      for( j = k+1,    iajk  = (k+1)+jak, jbj  = (k+1)*LDB;
           j < N; j++, iajk += 1,         jbj += LDB )
      {
         t0 = A[iajk];
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= t0 * B[ibik]; }
      }
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] *= ALPHA; }
   }
}

template<typename T>
static  void HPLAI_trsmRUNN
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   int64_t                        i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

   for( j = 0, jaj = 0, jbj = 0; j < N; j++, jaj += LDA, jbj += LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = 0, iakj = jaj, jbk = 0; k < j; k++, iakj += 1, jbk += LDB )
      {
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= A[iakj] * B[ibik]; }
      }
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] /= A[j+jaj]; }
   }
}

template<typename T>
static  void HPLAI_trsmRUNU
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   int64_t                        i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

   for( j = 0, jaj = 0, jbj = 0; j < N; j++, jaj += LDA, jbj += LDB )
   {
      for( i = 0, ibij = jbj; i < M; i++, ibij += 1 ) { B[ibij] *= ALPHA; }
      for( k = 0, iakj = jaj, jbk = 0; k < j; k++, iakj += 1, jbk += LDB )
      {
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= A[iakj] * B[ibik]; }
      }
   }
}

template<typename T>
static  void HPLAI_trsmRUTN
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   register T            t0;
   int64_t                        i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

   for( k = N-1,     jak  = (N-1)*LDA, jbk  = (N-1)*LDB;
        k >= 0; k--, jak -= LDA,       jbk -= LDB )
   {
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] /= A[k+jak]; }
      for( j = 0, iajk = jak, jbj = 0; j < k; j++, iajk += 1, jbj += LDB )
      {
         t0 = A[iajk];
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= t0 * B[ibik]; }
      }
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] *= ALPHA; }
   }
}

template<typename T>
static  void HPLAI_trsmRUTU
(
   const int64_t                  M,
   const int64_t                  N,
   const T               ALPHA,
   const T               * A,
   const int64_t                  LDA,
   T                     * B,
   const int64_t                  LDB
)
{
   register T            t0;
   int64_t                        i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

   for( k = N-1,     jak  = (N-1)*LDA, jbk  = (N-1)*LDB;
        k >= 0; k--, jak -= LDA,       jbk -= LDB )
   {
      for( j = 0, iajk = jak, jbj = 0; j < k; j++, iajk += 1, jbj += LDB )
      {
         t0 = A[iajk];
         for( i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1 )
         { B[ibij] -= t0 * B[ibik]; }
      }
      for( i = 0, ibik = jbk; i < M; i++, ibik += 1 ) { B[ibik] *= ALPHA; }
   }
}

template <>
void blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    blas::Side SIDE,
    blas::Uplo UPLO,
    blas::Op TRANS,
    blas::Diag DIAG,
    int64_t M,
    int64_t N,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> ALPHA,
    HPLAI_T_AFLOAT const *A,
    int64_t LDA,
    HPLAI_T_AFLOAT *B,
    int64_t LDB)
{
   //HPLAI_pabort( __LINE__, "blas::trsm", "Use HPLAI_GEN_BLASPP_TRSM" );
   if(layout != blas::Layout::ColMajor)
   {
      blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
         blas::Layout::ColMajor,
         ( SIDE == blas::Side::Right ? blas::Side::Left  : blas::Side::Right ),
         ( UPLO == blas::Uplo::Lower ? blas::Uplo::Upper : blas::Uplo::Lower ),
         TRANS, DIAG, N, M, ALPHA, A, LDA, B, LDB );
      return;
   }
   int64_t                        i, j;

   if( ( M == 0 ) || ( N == 0 ) ) return;
 
   if( ALPHA == HPLAI_rzero )
   {
      for( j = 0; j < N; j++ )
      {  for( i = 0; i < M; i++ ) *(B+i+j*LDB) = HPLAI_rzero; }
      return;
   }

   if( SIDE == blas::Side::Left )
   {
      if( UPLO == blas::Uplo::Upper )
      {
         if( TRANS == blas::Op::NoTrans )
         {
            if( DIAG == blas::Diag::NonUnit )
            {      HPLAI_trsmLUNN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPLAI_trsmLUNU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
         else
         {
            if( DIAG == blas::Diag::NonUnit )
            {      HPLAI_trsmLUTN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPLAI_trsmLUTU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
      }
      else
      {
         if( TRANS == blas::Op::NoTrans )
         {
            if( DIAG == blas::Diag::NonUnit )
            {      HPLAI_trsmLLNN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPLAI_trsmLLNU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
         else
         {
            if( DIAG == blas::Diag::NonUnit )
            {      HPLAI_trsmLLTN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPLAI_trsmLLTU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
      }
   }
   else
   {
      if( UPLO == blas::Uplo::Upper )
      {
         if( TRANS == blas::Op::NoTrans )
         {
            if( DIAG == blas::Diag::NonUnit )
            {      HPLAI_trsmRUNN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPLAI_trsmRUNU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
         else
         {
            if( DIAG == blas::Diag::NonUnit )
            {      HPLAI_trsmRUTN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPLAI_trsmRUTU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
      }
      else
      {
         if( TRANS == blas::Op::NoTrans )
         {
            if( DIAG == blas::Diag::NonUnit )
            {      HPLAI_trsmRLNN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPLAI_trsmRLNU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
         else
         {
            if( DIAG == blas::Diag::NonUnit )
            {      HPLAI_trsmRLTN( M, N, ALPHA, A, LDA, B, LDB ); }
            else { HPLAI_trsmRLTU( M, N, ALPHA, A, LDA, B, LDB ); }
         }
      }
   }
}

#endif

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

#ifndef HPLAI_GEN_BLASPP_TRSV

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


template<typename T>
static  void HPLAI_trsvLNN
(
   const int64_t                  N,
   const T               * A,
   const int64_t                  LDA,
   T                     * X,
   const int64_t                  INCX
)
{
   int64_t                        i, iaij, ix, j, jaj, jx, ldap1 = LDA + 1;
   register T            t0;

   for( j = 0, jaj = 0, jx  = 0; j < N; j++, jaj += ldap1, jx += INCX )
   {
      X[jx] /= A[jaj]; t0 = X[jx];
      for( i = j+1,    iaij  = jaj+1, ix  = jx + INCX;
           i < N; i++, iaij += 1,     ix += INCX ) { X[ix] -= t0 * A[iaij]; }
   }
}

template<typename T>
static  void HPLAI_trsvLNU
(
   const int64_t                  N,
   const T               * A,
   const int64_t                  LDA,
   T                     * X,
   const int64_t                  INCX
)
{
   int64_t                        i, iaij, ix, j, jaj, jx, ldap1 = LDA + 1;
   register T            t0;

   for( j = 0, jaj = 0, jx = 0; j < N; j++, jaj += ldap1, jx += INCX )
   {
      t0 = X[jx];
      for( i = j+1,    iaij  = jaj+1, ix  = jx + INCX;
           i < N; i++, iaij += 1,     ix += INCX ) { X[ix] -= t0 * A[iaij]; }
   }
}

template<typename T>
static  void HPLAI_trsvLTN
(
   const int64_t                  N,
   const T               * A,
   const int64_t                  LDA,
   T                     * X,
   const int64_t                  INCX
)
{
   int64_t                        i, iaij, ix, j, jaj, jx, ldap1 = LDA + 1;
   register T            t0;

   for( j = N-1,     jaj  = (N-1)*(ldap1), jx  = (N-1)*INCX;
        j >= 0; j--, jaj -= ldap1,         jx -= INCX )
   {
      t0 = X[jx];
      for( i = j+1,    iaij  = 1+jaj, ix  = jx + INCX;
           i < N; i++, iaij += 1,     ix += INCX ) { t0 -= A[iaij] * X[ix]; }
      t0 /= A[jaj]; X[jx] = t0;
   }
}

template<typename T>
static  void HPLAI_trsvLTU
(
   const int64_t                  N,
   const T               * A,
   const int64_t                  LDA,
   T                     * X,
   const int64_t                  INCX
)
{
   int64_t                        i, iaij, ix, j, jaj, jx, ldap1 = LDA + 1;
   register T            t0;

   for( j = N-1,     jaj  = (N-1)*(ldap1), jx  = (N-1)*INCX;
        j >= 0; j--, jaj -= ldap1,         jx -= INCX )
   {
      t0 = X[jx];
      for( i = j+1,    iaij  = 1+jaj, ix  = jx + INCX;
           i < N; i++, iaij += 1,     ix += INCX ) { t0 -= A[iaij] * X[ix]; }
      X[jx] = t0;
   }
}

template<typename T>
static  void HPLAI_trsvUNN
(
   const int64_t                  N,
   const T               * A,
   const int64_t                  LDA,
   T                     * X,
   const int64_t                  INCX
)
{
   int64_t                        i, iaij, ix, j, jaj, jx;
   register T            t0;

   for( j = N-1,     jaj  = (N-1)*LDA, jx  = (N-1)*INCX;
        j >= 0; j--, jaj -= LDA,       jx -= INCX )
   {
      X[jx] /= A[j+jaj]; t0 = X[jx];
      for( i = 0, iaij = jaj, ix = 0; i < j; i++, iaij += 1, ix += INCX )
      { X[ix] -= t0 * A[iaij]; }
   }
}

template<typename T>
static  void HPLAI_trsvUNU
(
   const int64_t                  N,
   const T               * A,
   const int64_t                  LDA,
   T                     * X,
   const int64_t                  INCX
)
{
   int64_t                        i, iaij, ix, j, jaj, jx;
   register T            t0;

   for( j = N-1,     jaj  = (N-1)*LDA, jx  = (N-1)*INCX;
        j >= 0; j--, jaj -= LDA,       jx -= INCX )
   {
      t0 = X[jx];
      for( i = 0, iaij = jaj, ix = 0; i < j; i++, iaij += 1, ix += INCX )
      { X[ix] -= t0 * A[iaij]; }
   }
}

template<typename T>
static  void HPLAI_trsvUTN
(
   const int64_t                  N,
   const T               * A,
   const int64_t                  LDA,
   T                     * X,
   const int64_t                  INCX
)
{
   int64_t                        i, iaij, ix, j, jaj, jx;
   register T            t0;

   for( j = 0, jaj = 0,jx = 0; j < N; j++, jaj += LDA, jx += INCX )
   {
      t0 = X[jx];
      for( i = 0, iaij = jaj, ix = 0; i < j; i++, iaij += 1, ix += INCX )
      { t0 -= A[iaij] * X[ix]; }
      t0 /= A[iaij]; X[jx] = t0;
   }
}

template<typename T>
static  void HPLAI_trsvUTU
(
   const int64_t                  N,
   const T               * A,
   const int64_t                  LDA,
   T                     * X,
   const int64_t                  INCX
)
{
   int64_t                        i, iaij, ix, j, jaj, jx;
   register T            t0;

   for( j = 0, jaj = 0, jx = 0; j < N; j++, jaj += LDA, jx += INCX )
   {
      t0 = X[jx];
      for( i = 0, iaij = jaj, ix = 0; i < j; i++, iaij += 1, ix += INCX )
      { t0 -= A[iaij] * X[ix]; }
      X[jx] = t0;
   }
}

template <>
void blas::trsv<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
	blas::Layout layout,
	blas::Uplo UPLO,
	blas::Op TRANS,
	blas::Diag DIAG,
	int64_t N,
	HPLAI_T_AFLOAT const *A,
	int64_t LDA,
	HPLAI_T_AFLOAT *X,
	int64_t INCX)
{
   //HPLAI_pabort( __LINE__, "blas::trsv", "Use Use HPLAI_GEN_BLASPP_TRSV" );
   if(layout != blas::Layout::ColMajor)
   {
      blas::trsv<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
         blas::Layout::ColMajor,
         ( UPLO  == blas::Uplo::Upper  ? blas::Uplo::Lower : blas::Uplo::Upper   ),
         ( TRANS == blas::Op::NoTrans ? blas::Op::Trans : blas::Op::NoTrans ),
         DIAG, N, A, LDA, X, INCX );
      return;
   }
   if( N == 0 ) return;
 
   if( UPLO == blas::Uplo::Upper )
   {
      if( TRANS == blas::Op::NoTrans )
      {
         if( DIAG == blas::Diag::NonUnit ) { HPLAI_trsvUNN( N,    A, LDA, X, INCX ); }
         else                     { HPLAI_trsvUNU( N,    A, LDA, X, INCX ); }
      }
      else
      {
         if( DIAG == blas::Diag::NonUnit ) { HPLAI_trsvUTN( N,    A, LDA, X, INCX ); }
         else                     { HPLAI_trsvUTU( N,    A, LDA, X, INCX ); }
      }
   }
   else
   {
      if( TRANS == blas::Op::NoTrans )
      {
         if( DIAG == blas::Diag::NonUnit ) { HPLAI_trsvLNN( N,    A, LDA, X, INCX ); }
         else                     { HPLAI_trsvLNU( N,    A, LDA, X, INCX ); }
      }
      else
      {
         if( DIAG == blas::Diag::NonUnit ) { HPLAI_trsvLTN( N,    A, LDA, X, INCX ); }
         else                     { HPLAI_trsvLTU( N,    A, LDA, X, INCX ); }
      }
   }
}

#endif