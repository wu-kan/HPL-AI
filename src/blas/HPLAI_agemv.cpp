/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef HPL_CALL_VSIPL

#ifdef STDC_HEADERS
static void HPLAI_agemv0
(
   const enum HPL_TRANS       TRANS,
   const int                  M,
   const int                  N,
   const HPLAI_T_AFLOAT               ALPHA,
   const HPLAI_T_AFLOAT               * A,
   const int                  LDA,
   const HPLAI_T_AFLOAT               * X,
   const int                  INCX,
   const HPLAI_T_AFLOAT               BETA,
   HPLAI_T_AFLOAT                     * Y,
   const int                  INCY
)
#else
static void HPLAI_agemv0( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
   const enum HPL_TRANS       TRANS;
   const int                  INCX, INCY, LDA, M, N;
   const HPLAI_T_AFLOAT               ALPHA, BETA;
   const HPLAI_T_AFLOAT               * A, * X;
   HPLAI_T_AFLOAT                     * Y;
#endif
{
/*
 * .. Local Variables ..
 */
   int                        i, iaij, ix, iy, j, jaj, jx, jy;
   register HPLAI_T_AFLOAT            t0;
/* ..
 * .. Executable Statements ..
 */
   if( ( M == 0 ) || ( N == 0 ) ||
       ( ( ALPHA == HPL_rzero ) && ( BETA == HPL_rone  ) ) ) return;
 
   if( ALPHA == HPL_rzero ) { HPLAI_ascal( M, BETA, Y, INCY ); return; }
 
   if( TRANS == HplNoTrans )
   {
      HPLAI_ascal( M, BETA, Y, INCY );
      for( j = 0, jaj  = 0, jx = 0; j < N; j++, jaj += LDA, jx += INCX )
      {
         t0 = ALPHA * X[jx];
         for( i = 0, iaij = jaj, iy = 0; i < M; i++, iaij += 1, iy += INCY )
         { Y[iy] += A[iaij] * t0; }
      }
   }
   else
   {
      for( j = 0, jaj  = 0, jy  = 0; j < N; j++, jaj += LDA, jy += INCY )
      {
         t0 = HPL_rzero;
         for( i = 0, iaij = jaj, ix = 0; i < M; i++, iaij += 1, ix += INCX )
         { t0 += A[iaij] * X[ix]; }
         if( BETA == HPL_rzero ) Y[jy] = ALPHA * t0;
         else                    Y[jy] = BETA * Y[jy] + ALPHA * t0;
      }
   }
}
#endif

#ifdef STDC_HEADERS
void HPLAI_agemv
(
   const enum HPL_ORDER             ORDER,
   const enum HPL_TRANS             TRANS,
   const int                        M,
   const int                        N,
   const HPLAI_T_AFLOAT                     ALPHA,
   const HPLAI_T_AFLOAT *                   A,
   const int                        LDA,
   const HPLAI_T_AFLOAT *                   X,
   const int                        INCX,
   const HPLAI_T_AFLOAT                     BETA,
   HPLAI_T_AFLOAT *                         Y,
   const int                        INCY
)
#else
void HPLAI_agemv
( ORDER, TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
   const enum HPL_ORDER             ORDER;
   const enum HPL_TRANS             TRANS;
   const int                        M;
   const int                        N;
   const HPLAI_T_AFLOAT                     ALPHA;
   const HPLAI_T_AFLOAT *                   A;
   const int                        LDA;
   const HPLAI_T_AFLOAT *                   X;
   const int                        INCX;
   const HPLAI_T_AFLOAT                     BETA;
   HPLAI_T_AFLOAT *                         Y;
   const int                        INCY;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_agemv performs one of the matrix-vector operations
 *  
 *     y := alpha * op( A ) * x + beta * y,
 *  
 *  where op( X ) is one of
 *  
 *     op( X ) = X   or   op( X ) = X^T.
 *  
 * where alpha and beta are scalars, x and y are vectors and  A  is an m
 * by n matrix.
 *
 * Arguments
 * =========
 *
 * ORDER   (local input)                 const enum HPL_ORDER
 *         On entry, ORDER  specifies the storage format of the operands
 *         as follows:                                                  
 *            ORDER = HplRowMajor,                                      
 *            ORDER = HplColumnMajor.                                   
 *
 * TRANS   (local input)                 const enum HPL_TRANS
 *         On entry,  TRANS  specifies the  operation to be performed as
 *         follows:   
 *            TRANS = HplNoTrans y := alpha*A  *x + beta*y,
 *            TRANS = HplTrans   y := alpha*A^T*x + beta*y.
 *
 * M       (local input)                 const int
 *         On entry,  M  specifies  the number of rows of  the matrix A.
 *         M must be at least zero.
 *
 * N       (local input)                 const int
 *         On entry, N  specifies the number of columns of the matrix A.
 *         N must be at least zero.
 *
 * ALPHA   (local input)                 const HPLAI_T_AFLOAT
 *         On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *         supplied as zero then  A and X  need not be set on input.
 *
 * A       (local input)                 const HPLAI_T_AFLOAT *
 *         On entry,  A  points  to an array of size equal to or greater
 *         than LDA * n.  Before  entry, the leading m by n part  of the
 *         array  A  must contain the matrix coefficients.
 *
 * LDA     (local input)                 const int
 *         On entry,  LDA  specifies  the  leading  dimension  of  A  as
 *         declared  in  the  calling  (sub) program.  LDA  must  be  at
 *         least MAX(1,m).
 *
 * X       (local input)                 const HPLAI_T_AFLOAT *
 *         On entry,  X  is an incremented array of dimension  at  least
 *         ( 1 + ( n - 1 ) * abs( INCX ) )  that  contains the vector x.
 *
 * INCX    (local input)                 const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero.
 *
 * BETA    (local input)                 const HPLAI_T_AFLOAT
 *         On entry, BETA  specifies the scalar beta.    When  ALPHA  is
 *         supplied as zero then  Y  need not be set on input.
 *
 * Y       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry,  Y  is an incremented array of dimension  at  least
 *         ( 1 + ( n - 1 ) * abs( INCY ) )  that  contains the vector y.
 *         Before entry with BETA non-zero, the incremented array Y must
 *         contain the vector  y.  On exit,  Y  is  overwritten  by  the
 *         updated vector y.
 *
 * INCY    (local input)                 const int
 *         On entry, INCY specifies the increment for the elements of Y.
 *         INCY must not be zero.
 *
 * ---------------------------------------------------------------------
 */ 
#ifdef HPL_CALL_CBLAS
   cblas_sgemv( ORDER, TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY );
#endif
#ifdef HPL_CALL_VSIPL
   if( ORDER == HplColumnMajor )
   {
      HPLAI_agemv0( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY );
   }
   else
   {
      HPLAI_agemv0( ( TRANS == HplNoTrans ? HplTrans : HplNoTrans ),
                  N, M, ALPHA, A, LDA, X, INCX, BETA, Y, INCY );
   }
#endif
#ifdef HPL_CALL_FBLAS
   HPLAI_T_AFLOAT                    alpha = ALPHA, beta = BETA;
#ifdef StringSunStyle
#ifdef HPL_USE_F77_INTEGER_DEF
   F77_INTEGER               IONE = 1;
#else
   int                       IONE = 1;
#endif
#endif
#ifdef StringStructVal
   F77_CHAR                  ftran;
#endif
#ifdef StringStructPtr
   F77_CHAR                  ftran;
#endif
#ifdef StringCrayStyle
   F77_CHAR                  ftran;
#endif
 
#ifdef HPL_USE_F77_INTEGER_DEF
   const F77_INTEGER         F77M    = M,   F77N    = N,
                             F77lda  = LDA, F77incx = INCX, F77incy = INCY;
#else
#define F77M                 M
#define F77N                 N
#define F77lda               LDA
#define F77incx              INCX
#define F77incy              INCY
#endif
   char                      ctran;

   if( ORDER == HplColumnMajor )
   {
      ctran = ( TRANS == HplNoTrans ? 'N' : 'T' );

#ifdef StringSunStyle
      F77dgemv( &ctran, &F77M, &F77N, &alpha, A, &F77lda, X, &F77incx,
                &beta, Y, &F77incy, IONE );
#endif
#ifdef StringCrayStyle
      ftran = HPL_C2F_CHAR( ctran );
      F77dgemv( ftran,  &F77M, &F77N, &alpha, A, &F77lda, X, &F77incx,
                &beta, Y, &F77incy );
#endif
#ifdef StringStructVal
      ftran.len = 1; ftran.cp = &ctran;
      F77dgemv( ftran,  &F77M, &F77N, &alpha, A, &F77lda, X, &F77incx,
                &beta, Y, &F77incy );
#endif
#ifdef StringStructPtr
      ftran.len = 1; ftran.cp = &ctran;
      F77dgemv( &ftran, &F77M, &F77N, &alpha, A, &F77lda, X, &F77incx,
                &beta, Y, &F77incy );
#endif
   }
   else
   {
      ctran = ( TRANS == HplNoTrans ? 'T' : 'N' );
#ifdef StringSunStyle
      F77dgemv( &ctran, &F77N, &F77M, &alpha, A, &F77lda, X, &F77incx,
                &beta, Y, &F77incy, IONE );
#endif
#ifdef StringCrayStyle
      ftran = HPL_C2F_CHAR( ctran );
      F77dgemv( ftran,  &F77N, &F77M, &alpha, A, &F77lda, X, &F77incx,
                &beta, Y, &F77incy );
#endif
#ifdef StringStructVal
      ftran.len = 1; ftran.cp = &ctran;
      F77dgemv( ftran,  &F77N, &F77M, &alpha, A, &F77lda, X, &F77incx,
                &beta, Y, &F77incy );
#endif
#ifdef StringStructPtr
      ftran.len = 1; ftran.cp = &ctran;
      F77dgemv( &ftran, &F77N, &F77M, &alpha, A, &F77lda, X, &F77incx,
                &beta, Y, &F77incy );
#endif
   }

#endif
/*
 * End of HPLAI_agemv
 */
}

#ifdef __cplusplus
}
#endif
