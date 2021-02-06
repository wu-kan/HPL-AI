/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPLAI_ager
(
   const enum HPL_ORDER             ORDER,
   const int                        M,
   const int                        N,
   const HPLAI_T_AFLOAT                     ALPHA,
   const HPLAI_T_AFLOAT *                   X,
   const int                        INCX,
   HPLAI_T_AFLOAT *                         Y,
   const int                        INCY,
   HPLAI_T_AFLOAT *                         A,
   const int                        LDA
)
#else
void HPLAI_ager
( ORDER, M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
   const enum HPL_ORDER             ORDER;
   const int                        M;
   const int                        N;
   const HPLAI_T_AFLOAT                     ALPHA;
   const HPLAI_T_AFLOAT *                   X;
   const int                        INCX;
   HPLAI_T_AFLOAT *                         Y;
   const int                        INCY;
   HPLAI_T_AFLOAT *                         A;
   const int                        LDA;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_ager performs the rank 1 operation
 *  
 *     A := alpha * x * y^T + A,
 *  
 * where alpha is a scalar,  x is an m-element vector, y is an n-element
 * vector and A is an m by n matrix.
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
 *         supplied as zero then  X and Y  need not be set on input.
 *
 * X       (local input)                 const HPLAI_T_AFLOAT *
 *         On entry,  X  is an incremented array of dimension  at  least
 *         ( 1 + ( m - 1 ) * abs( INCX ) )  that  contains the vector x.
 *
 * INCX    (local input)                 const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero.
 *
 * Y       (local input)                 HPLAI_T_AFLOAT *
 *         On entry,  Y  is an incremented array of dimension  at  least
 *         ( 1 + ( n - 1 ) * abs( INCY ) )  that  contains the vector y.
 *
 * INCY    (local input)                 const int
 *         On entry, INCY specifies the increment for the elements of Y.
 *         INCY must not be zero.
 *
 * A       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry,  A  points  to an array of size equal to or greater
 *         than LDA * n.  Before  entry, the leading m by n part  of the
 *         array  A  must contain the matrix coefficients. On exit, A is
 *         overwritten by the updated matrix.
 *
 * LDA     (local input)                 const int
 *         On entry,  LDA  specifies  the  leading  dimension  of  A  as
 *         declared  in  the  calling  (sub) program.  LDA  must  be  at
 *         least MAX(1,m).
 *
 * ---------------------------------------------------------------------
 */ 
#ifdef HPL_CALL_CBLAS
   cblas_sger( ORDER, M, N, ALPHA, X, INCX, Y, INCY, A, LDA );
#endif
#ifdef HPL_CALL_VSIPL
   register HPLAI_T_AFLOAT           t0;
   int                       i, iaij, ix, iy, j, jaj, jx, jy;

   if( ( M == 0 ) || ( N == 0 ) || ( ALPHA == HPL_rzero ) ) return;
 
   if( ORDER == HplColumnMajor )
   {
      for( j = 0, jaj = 0, jy = 0; j < N; j++, jaj += LDA, jy += INCY )
      {
         t0 = ALPHA * Y[jy];
         for( i = 0, iaij = jaj, ix = 0; i < M; i++, iaij += 1, ix += INCX )
         { A[iaij] += X[ix] * t0; }
      }
   }
   else
   {
      for( j = 0, jaj = 0, jx = 0; j < M; j++, jaj += LDA, jx += INCX )
      {
         t0 = ALPHA * X[jx];
         for( i = 0, iaij = jaj, iy = 0; i < N; i++, iaij += 1, iy += INCY )
         { A[iaij] += Y[iy] * t0; }
      }
   }
#endif
#ifdef HPL_CALL_FBLAS
   HPLAI_T_AFLOAT                    alpha = ALPHA;
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

   if( ORDER == HplColumnMajor )
   {  F77dger( &F77M, &F77N, &alpha, X, &F77incx, Y, &F77incy, A, &F77lda ); }
   else
   {  F77dger( &F77N, &F77M, &alpha, Y, &F77incy, X, &F77incx, A, &F77lda ); }
#endif
/*
 * End of HPLAI_ager
 */
}

#ifdef __cplusplus
}
#endif
