/*
 * Include files
 */
#include "hpl_ai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef HPL_CALL_VSIPL

#ifdef STDC_HEADERS
static void HPL_AI_agemmNN
(
   const int                  M,
   const int                  N,
   const int                  K,
   const HPL_AI_T_afloat               ALPHA,
   const HPL_AI_T_afloat               * A,
   const int                  LDA,
   const HPL_AI_T_afloat               * B,
   const int                  LDB,
   const HPL_AI_T_afloat               BETA,
   HPL_AI_T_afloat                     * C,
   const int                  LDC
)
#else
static void HPL_AI_agemmNN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
   const int                  K, LDA, LDB, LDC, M, N;
   const HPL_AI_T_afloat               ALPHA, BETA;
   const HPL_AI_T_afloat               * A, * B;
   HPL_AI_T_afloat                     * C;
#endif
{
   register HPL_AI_T_afloat            t0;
   int                        i, iail, iblj, icij, j, jal, jbj, jcj, l;

   for( j = 0, jbj = 0, jcj  = 0; j < N; j++, jbj += LDB, jcj += LDC )
   {
      HPL_AI_ascal( M, BETA, C+jcj, 1 );
      for( l = 0, jal = 0, iblj = jbj; l < K; l++, jal += LDA, iblj += 1 )
      {
         t0 = ALPHA * B[iblj];
         for( i = 0, iail = jal, icij = jcj; i < M; i++, iail += 1, icij += 1 )
         { C[icij] += A[iail] * t0; }
      }
   }
}

#ifdef STDC_HEADERS
static void HPL_AI_agemmNT
(
   const int                  M,
   const int                  N,
   const int                  K,
   const HPL_AI_T_afloat               ALPHA,
   const HPL_AI_T_afloat               * A,
   const int                  LDA,
   const HPL_AI_T_afloat               * B,
   const int                  LDB,
   const HPL_AI_T_afloat               BETA,
   HPL_AI_T_afloat                     * C,
   const int                  LDC
)
#else
static void HPL_AI_agemmNT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
   const int                  K, LDA, LDB, LDC, M, N;
   const HPL_AI_T_afloat               ALPHA, BETA;
   const HPL_AI_T_afloat               * A, * B;
   HPL_AI_T_afloat                     * C;
#endif
{
   register HPL_AI_T_afloat            t0;
   int                        i, iail, ibj, ibjl, icij, j, jal, jcj, l;

   for( j = 0, ibj  = 0, jcj  = 0; j < N; j++, ibj += 1, jcj += LDC )
   {
      HPL_AI_ascal( M, BETA, C+jcj, 1 );
      for( l = 0, jal = 0, ibjl = ibj; l < K; l++, jal += LDA, ibjl += LDB )
      {
         t0 = ALPHA * B[ibjl];
         for( i = 0, iail = jal, icij = jcj; i < M; i++, iail += 1, icij += 1 )
         { C[icij] += A[iail] * t0; }
      }
   }
}

#ifdef STDC_HEADERS
static void HPL_AI_agemmTN
(
   const int                  M,
   const int                  N,
   const int                  K,
   const HPL_AI_T_afloat               ALPHA,
   const HPL_AI_T_afloat               * A,
   const int                  LDA,
   const HPL_AI_T_afloat               * B,
   const int                  LDB,
   const HPL_AI_T_afloat               BETA,
   HPL_AI_T_afloat                     * C,
   const int                  LDC
)
#else
static void HPL_AI_agemmTN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
   const int                  K, LDA, LDB, LDC, M, N;
   const HPL_AI_T_afloat               ALPHA, BETA;
   const HPL_AI_T_afloat               * A, * B;
   HPL_AI_T_afloat                     * C;
#endif
{
   register HPL_AI_T_afloat            t0;
   int                        i, iai, iail, iblj, icij, j, jbj, jcj, l;

   for( j = 0, jbj = 0, jcj = 0; j < N; j++, jbj += LDB, jcj += LDC )
   {
      for( i = 0, icij = jcj, iai = 0; i < M; i++, icij += 1, iai += LDA )
      {
         t0 = HPL_rzero;
         for( l = 0, iail = iai, iblj = jbj; l < K; l++, iail += 1, iblj += 1 )
         { t0 += A[iail] * B[iblj]; }
         if( BETA == HPL_rzero ) C[icij]  = HPL_rzero;
         else                    C[icij] *= BETA;
         C[icij] += ALPHA * t0;
      }
   }
}

#ifdef STDC_HEADERS
static void HPL_AI_agemmTT
(
   const int                  M,
   const int                  N,
   const int                  K,
   const HPL_AI_T_afloat               ALPHA,
   const HPL_AI_T_afloat               * A,
   const int                  LDA,
   const HPL_AI_T_afloat               * B,
   const int                  LDB,
   const HPL_AI_T_afloat               BETA,
   HPL_AI_T_afloat                     * C,
   const int                  LDC
)
#else
static void HPL_AI_agemmTT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
   const int                  K, LDA, LDB, LDC, M, N;
   const HPL_AI_T_afloat               ALPHA, BETA;
   const HPL_AI_T_afloat               * A, * B;
   HPL_AI_T_afloat                     * C;
#endif
{
   register HPL_AI_T_afloat            t0;
   int                        i, iali, ibj, ibjl, icij, j, jai, jcj, l;

   for( j = 0, ibj = 0, jcj  = 0; j < N; j++, ibj += 1, jcj += LDC )
   {
      for( i = 0, icij = jcj, jai = 0; i < M; i++, icij += 1, jai += LDA )
      {
         t0 = HPL_rzero;
         for( l = 0,      iali  = jai, ibjl  = ibj;
              l < K; l++, iali += 1,   ibjl += LDB ) t0 += A[iali] * B[ibjl];
         if( BETA == HPL_rzero ) C[icij]  = HPL_rzero;
         else                    C[icij] *= BETA;
         C[icij] += ALPHA * t0;
      }
   }
}

#ifdef STDC_HEADERS
static void HPL_AI_agemm0
(
   const enum HPL_TRANS       TRANSA,
   const enum HPL_TRANS       TRANSB,
   const int                  M,
   const int                  N,
   const int                  K,
   const HPL_AI_T_afloat               ALPHA,
   const HPL_AI_T_afloat               * A,
   const int                  LDA,
   const HPL_AI_T_afloat               * B,
   const int                  LDB,
   const HPL_AI_T_afloat               BETA,
   HPL_AI_T_afloat                     * C,
   const int                  LDC
)
#else
static void HPL_AI_agemm0( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
                        BETA, C, LDC )
   const enum HPL_TRANS       TRANSA, TRANSB;
   const int                  K, LDA, LDB, LDC, M, N;
   const HPL_AI_T_afloat               ALPHA, BETA;
   const HPL_AI_T_afloat               * A, * B;
   HPL_AI_T_afloat                     * C;
#endif
{
   int                        i, j;

   if( ( M == 0 ) || ( N == 0 ) ||
       ( ( ( ALPHA == HPL_rzero ) || ( K == 0 ) ) &&
         ( BETA == HPL_rone ) ) ) return;

   if( ALPHA == HPL_rzero )
   {
      for( j = 0; j < N; j++ )
      {  for( i = 0; i < M; i++ ) *(C+i+j*LDC) = HPL_rzero; }
      return;
   }

   if( TRANSB == HplNoTrans )
   {
      if( TRANSA == HplNoTrans )
      { HPL_AI_agemmNN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
      else
      { HPL_AI_agemmTN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
   }
   else
   {
      if( TRANSA == HplNoTrans )
      { HPL_AI_agemmNT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
      else
      { HPL_AI_agemmTT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ); }
   }
}

#endif

#ifdef STDC_HEADERS
void HPL_AI_agemm
(
   const enum HPL_ORDER             ORDER,
   const enum HPL_TRANS             TRANSA,
   const enum HPL_TRANS             TRANSB,
   const int                        M,
   const int                        N,
   const int                        K,
   const HPL_AI_T_afloat                     ALPHA,
   const HPL_AI_T_afloat *                   A,
   const int                        LDA,
   const HPL_AI_T_afloat *                   B,
   const int                        LDB,
   const HPL_AI_T_afloat                     BETA,
   HPL_AI_T_afloat *                         C,
   const int                        LDC
)
#else
void HPL_AI_agemm
( ORDER, TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
   const enum HPL_ORDER             ORDER;
   const enum HPL_TRANS             TRANSA;
   const enum HPL_TRANS             TRANSB;
   const int                        M;
   const int                        N;
   const int                        K;
   const HPL_AI_T_afloat                     ALPHA;
   const HPL_AI_T_afloat *                   A;
   const int                        LDA;
   const HPL_AI_T_afloat *                   B;
   const int                        LDB;
   const HPL_AI_T_afloat                     BETA;
   HPL_AI_T_afloat *                         C;
   const int                        LDC;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_AI_agemm performs one of the matrix-matrix operations
 *  
 *     C := alpha * op( A ) * op( B ) + beta * C
 *  
 *  where op( X ) is one of
 *  
 *     op( X ) = X   or   op( X ) = X^T.
 *  
 * Alpha and beta are scalars,  and A,  B and C are matrices, with op(A)
 * an m by k matrix, op(B) a k by n matrix and  C an m by n matrix.
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
 * TRANSA  (local input)                 const enum HPL_TRANS
 *         On entry, TRANSA  specifies the form of  op(A)  to be used in
 *         the matrix-matrix operation follows:                         
 *            TRANSA==HplNoTrans    : op( A ) = A,                     
 *            TRANSA==HplTrans      : op( A ) = A^T,                   
 *            TRANSA==HplConjTrans  : op( A ) = A^T.                   
 *
 * TRANSB  (local input)                 const enum HPL_TRANS
 *         On entry, TRANSB  specifies the form of  op(B)  to be used in
 *         the matrix-matrix operation follows:                         
 *            TRANSB==HplNoTrans    : op( B ) = B,                     
 *            TRANSB==HplTrans      : op( B ) = B^T,                   
 *            TRANSB==HplConjTrans  : op( B ) = B^T.                   
 *
 * M       (local input)                 const int
 *         On entry,  M  specifies  the  number  of rows  of the  matrix
 *         op(A)  and  of  the  matrix  C.  M  must  be  at least  zero.
 *
 * N       (local input)                 const int
 *         On entry,  N  specifies  the number  of columns of the matrix
 *         op(B)  and  the number of columns of the matrix  C. N must be
 *         at least zero.
 *
 * K       (local input)                 const int
 *         On entry,  K  specifies  the  number of columns of the matrix
 *         op(A) and the number of rows of the matrix op(B).  K  must be
 *         be at least  zero.
 *
 * ALPHA   (local input)                 const HPL_AI_T_afloat
 *         On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *         supplied  as  zero  then the elements of the matrices A and B
 *         need not be set on input.
 *
 * A       (local input)                 const HPL_AI_T_afloat *
 *         On entry,  A  is an array of dimension (LDA,ka),  where ka is
 *         k  when   TRANSA==HplNoTrans,  and  is  m  otherwise.  Before
 *         entry  with  TRANSA==HplNoTrans, the  leading  m by k part of
 *         the array  A must contain the matrix A, otherwise the leading
 *         k  by  m  part of the array  A  must  contain the  matrix  A.
 *
 * LDA     (local input)                 const int
 *         On entry, LDA  specifies the first dimension of A as declared
 *         in the  calling (sub) program. When  TRANSA==HplNoTrans  then
 *         LDA must be at least max(1,m), otherwise LDA must be at least
 *         max(1,k).
 *
 * B       (local input)                 const HPL_AI_T_afloat *
 *         On entry, B is an array of dimension (LDB,kb),  where  kb  is
 *         n   when  TRANSB==HplNoTrans, and  is  k  otherwise.   Before
 *         entry with TRANSB==HplNoTrans,  the  leading  k by n  part of
 *         the array  B must contain the matrix B, otherwise the leading
 *         n  by  k  part of the array  B  must  contain  the matrix  B.
 *
 * LDB     (local input)                 const int
 *         On entry, LDB  specifies the first dimension of B as declared
 *         in the  calling (sub) program. When  TRANSB==HplNoTrans  then
 *         LDB must be at least max(1,k), otherwise LDB must be at least
 *         max(1,n).
 *
 * BETA    (local input)                 const HPL_AI_T_afloat
 *         On entry,  BETA  specifies the scalar  beta.   When  BETA  is
 *         supplied  as  zero  then  the  elements of the matrix C  need
 *         not be set on input.
 *
 * C       (local input/output)          HPL_AI_T_afloat *
 *         On entry,  C  is an array of dimension (LDC,n). Before entry,
 *         the  leading m by n part  of  the  array  C  must contain the
 *         matrix C,  except when beta is zero, in which case C need not
 *         be set on entry. On exit, the array  C  is overwritten by the
 *         m by n  matrix ( alpha*op( A )*op( B ) + beta*C ).
 *
 * LDC     (local input)                 const int
 *         On entry, LDC  specifies the first dimension of C as declared
 *         in  the   calling  (sub)  program.   LDC  must  be  at  least
 *         max(1,m).
 *
 * ---------------------------------------------------------------------
 */ 
#ifdef HPL_CALL_CBLAS
   cblas_sgemm( ORDER, TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
                BETA, C, LDC );
#endif
#ifdef HPL_CALL_VSIPL
   if( ORDER == HplColumnMajor )
   {
      HPL_AI_agemm0( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA,
                  C, LDC );
   }
   else
   {
      HPL_AI_agemm0( TRANSB, TRANSA, N, M, K, ALPHA, B, LDB, A, LDA, BETA,
                  C, LDC );
   }
#endif
#ifdef HPL_CALL_FBLAS
   HPL_AI_T_afloat                    alpha = ALPHA, beta = BETA;
#ifdef StringSunStyle
#ifdef HPL_USE_F77_INTEGER_DEF
   F77_INTEGER               IONE = 1;
#else
   int                       IONE = 1;
#endif
#endif
#ifdef StringStructVal
   F77_CHAR                  ftransa;
   F77_CHAR                  ftransb;
#endif
#ifdef StringStructPtr
   F77_CHAR                  ftransa;
   F77_CHAR                  ftransb;
#endif
#ifdef StringCrayStyle
   F77_CHAR                  ftransa;
   F77_CHAR                  ftransb;
#endif
#ifdef HPL_USE_F77_INTEGER_DEF
   const F77_INTEGER         F77M   = M,   F77N   = N,   F77K = K,
                             F77lda = LDA, F77ldb = LDB, F77ldc = LDC;
#else
#define F77M                 M
#define F77N                 N
#define F77K                 K
#define F77lda               LDA
#define F77ldb               LDB
#define F77ldc               LDC
#endif
   char                      ctransa, ctransb;

   if(      TRANSA == HplNoTrans ) ctransa = 'N';
   else if( TRANSA == HplTrans   ) ctransa = 'T';
   else                            ctransa = 'C';
 
   if(      TRANSB == HplNoTrans ) ctransb = 'N';
   else if( TRANSB == HplTrans   ) ctransb = 'T';
   else                            ctransb = 'C';

   if( ORDER == HplColumnMajor )
   {
#ifdef StringSunStyle
      F77dgemm( &ctransa, &ctransb, &F77M, &F77N, &F77K, &alpha, A, &F77lda,
                B, &F77ldb, &beta, C, &F77ldc, IONE, IONE );
#endif
#ifdef StringCrayStyle
      ftransa = HPL_C2F_CHAR( ctransa ); ftransb = HPL_C2F_CHAR( ctransb );
      F77dgemm( ftransa,  ftransb,  &F77M, &F77N, &F77K, &alpha, A, &F77lda,
                B, &F77ldb, &beta, C, &F77ldc );
#endif
#ifdef StringStructVal
      ftransa.len = 1; ftransa.cp = &ctransa;
      ftransb.len = 1; ftransb.cp = &ctransb;
      F77dgemm( ftransa,  ftransb,  &F77M, &F77N, &F77K, &alpha, A, &F77lda,
                B, &F77ldb, &beta, C, &F77ldc );
#endif
#ifdef StringStructPtr
      ftransa.len = 1; ftransa.cp = &ctransa;
      ftransb.len = 1; ftransb.cp = &ctransb;
      F77dgemm( &ftransa, &ftransb, &F77M, &F77N, &F77K, &alpha, A, &F77lda,
                B, &F77ldb, &beta, C, &F77ldc );
#endif
   }
   else
   {
#ifdef StringSunStyle
      F77dgemm( &ctransb, &ctransa, &F77N, &F77M, &F77K, &alpha, B, &F77ldb,
                A, &F77lda, &beta, C, &F77ldc, IONE, IONE );
#endif
#ifdef StringCrayStyle
      ftransa = HPL_C2F_CHAR( ctransa ); ftransb = HPL_C2F_CHAR( ctransb );
      F77dgemm( ftransb,  ftransa,  &F77N, &F77M, &F77K, &alpha, B, &F77ldb,
                A, &F77lda, &beta, C, &F77ldc );
#endif
#ifdef StringStructVal
      ftransa.len = 1; ftransa.cp = &ctransa;
      ftransb.len = 1; ftransb.cp = &ctransb;
      F77dgemm( ftransb,  ftransa,  &F77N, &F77M, &F77K, &alpha, B, &F77ldb,
                A, &F77lda, &beta, C, &F77ldc );
#endif
#ifdef StringStructPtr
      ftransa.len = 1; ftransa.cp = &ctransa;
      ftransb.len = 1; ftransb.cp = &ctransb;
      F77dgemm( &ftransb, &ftransa, &F77N, &F77M, &F77K, &alpha, B, &F77ldb,
                A, &F77lda, &beta, C, &F77ldc );
#endif
   }
#endif
/*
 * End of HPL_AI_agemm
 */
}

#ifdef __cplusplus
}
#endif
