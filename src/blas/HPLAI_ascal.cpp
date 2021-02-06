/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPLAI_ascal
(
   const int                        N,
   const HPLAI_T_AFLOAT                     ALPHA,
   HPLAI_T_AFLOAT *                         X,
   const int                        INCX
)
#else
void HPLAI_ascal
( N, ALPHA, X, INCX )
   const int                        N;
   const HPLAI_T_AFLOAT                     ALPHA;
   HPLAI_T_AFLOAT *                         X;
   const int                        INCX;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_ascal scales the vector x by alpha.
 * 
 *
 * Arguments
 * =========
 *
 * N       (local input)                 const int
 *         On entry, N specifies the length of the vector x. N  must  be
 *         at least zero.
 *
 * ALPHA   (local input)                 const HPLAI_T_AFLOAT
 *         On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *         supplied as zero, then the entries of the incremented array X
 *         need not be set on input.
 *
 * X       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry,  X  is an incremented array of dimension  at  least
 *         ( 1 + ( n - 1 ) * abs( INCX ) )  that  contains the vector x.
 *         On exit, the entries of the incremented array  X  are  scaled
 *         by the scalar alpha.
 *
 * INCX    (local input)                 const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero.
 *
 * ---------------------------------------------------------------------
 */ 
#ifdef HPL_CALL_CBLAS
   cblas_sscal( N, ALPHA, X, INCX );
#endif
#ifdef HPL_CALL_VSIPL
   register HPLAI_T_AFLOAT           x0, x1, x2, x3, x4, x5, x6, x7;
   register const HPLAI_T_AFLOAT     alpha = ALPHA;
   const HPLAI_T_AFLOAT              * StX;
   register int              i;
   int                       nu;
   const int                 incX2 = 2 * INCX, incX3 = 3 * INCX,
                             incX4 = 4 * INCX, incX5 = 5 * INCX,
                             incX6 = 6 * INCX, incX7 = 7 * INCX,
                             incX8 = 8 * INCX;

   if( ( N > 0 ) && ( alpha != HPL_rone ) )
   {
      if( alpha == HPL_rzero )
      {
         if( ( nu = ( N >> 3 ) << 3 ) != 0 )
         {
            StX = (HPLAI_T_AFLOAT *)X + nu * INCX;
 
            do
            {
               (*X)     = HPL_rzero; X[incX4] = HPL_rzero;
               X[INCX ] = HPL_rzero; X[incX5] = HPL_rzero;
               X[incX2] = HPL_rzero; X[incX6] = HPL_rzero;
               X[incX3] = HPL_rzero; X[incX7] = HPL_rzero; X += incX8;

            } while( X != StX );
         }
 
         for( i = N - nu; i != 0; i-- ) { *X = HPL_rzero; X += INCX; }
      }
      else
      {
         if( ( nu = ( N >> 3 ) << 3 ) != 0 )
         {
            StX = X + nu * INCX;
 
            do
            {
               x0 = (*X);     x4 = X[incX4]; x1 = X[INCX ]; x5 = X[incX5];
               x2 = X[incX2]; x6 = X[incX6]; x3 = X[incX3]; x7 = X[incX7];
 
               x0 *= alpha;   x4 *= alpha;   x1 *= alpha;   x5 *= alpha;
               x2 *= alpha;   x6 *= alpha;   x3 *= alpha;   x7 *= alpha;
 
               (*X)     = x0; X[incX4] = x4; X[INCX ] = x1; X[incX5] = x5;
               X[incX2] = x2; X[incX6] = x6; X[incX3] = x3; X[incX7] = x7;
 
               X  += incX8;
 
            } while( X != StX );
         }
 
         for( i = N - nu; i != 0; i-- )
         { x0 = (*X); x0 *= alpha; *X = x0; X += INCX; }
      }
   }
#endif
#ifdef HPL_CALL_FBLAS
   HPLAI_T_AFLOAT                    alpha = ALPHA;
#ifdef HPL_USE_F77_INTEGER_DEF
   const F77_INTEGER         F77N = N, F77incx = INCX;
#else
#define F77N                 N
#define F77incx              INCX
#endif

   F77dscal( &F77N, &alpha, X, &F77incx );
#endif
/*
 * End of HPLAI_ascal
 */
}
 
#ifdef __cplusplus
}
#endif
