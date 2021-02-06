/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPLAI_acopy
(
   const int                        N,
   const HPLAI_T_AFLOAT *                   X,
   const int                        INCX,
   HPLAI_T_AFLOAT *                         Y,
   const int                        INCY
)
#else
void HPLAI_acopy
( N, X, INCX, Y, INCY )
   const int                        N;
   const HPLAI_T_AFLOAT *                   X;
   const int                        INCX;
   HPLAI_T_AFLOAT *                         Y;
   const int                        INCY;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_acopy copies the vector x into the vector y.
 * 
 *
 * Arguments
 * =========
 *
 * N       (local input)                 const int
 *         On entry, N specifies the length of the vectors  x  and  y. N
 *         must be at least zero.
 *
 * X       (local input)                 const HPLAI_T_AFLOAT *
 *         On entry,  X  is an incremented array of dimension  at  least
 *         ( 1 + ( n - 1 ) * abs( INCX ) )  that  contains the vector x.
 *
 * INCX    (local input)                 const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero.
 *
 * Y       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry,  Y  is an incremented array of dimension  at  least
 *         ( 1 + ( n - 1 ) * abs( INCY ) )  that  contains the vector y.
 *         On exit, the entries of the incremented array  Y  are updated
 *         with the entries of the incremented array X.
 *
 * INCY    (local input)                 const int
 *         On entry, INCY specifies the increment for the elements of Y.
 *         INCY must not be zero.
 *
 * ---------------------------------------------------------------------
 */ 
#ifdef HPL_CALL_CBLAS
   cblas_scopy( N, X, INCX, Y, INCY );
#endif
#ifdef HPL_CALL_VSIPL
   register HPLAI_T_AFLOAT           x0, x1, x2, x3, x4, x5, x6, x7;
   const HPLAI_T_AFLOAT              * StX;
   register int              i;
   int                       nu;
   const int                 incX2 = 2 * INCX, incY2 = 2 * INCY,
                             incX3 = 3 * INCX, incY3 = 3 * INCY,
                             incX4 = 4 * INCX, incY4 = 4 * INCY,
                             incX5 = 5 * INCX, incY5 = 5 * INCY,
                             incX6 = 6 * INCX, incY6 = 6 * INCY,
                             incX7 = 7 * INCX, incY7 = 7 * INCY,
                             incX8 = 8 * INCX, incY8 = 8 * INCY;

   if( N > 0 )
   {
      if( ( nu = ( N >> 3 ) << 3 ) != 0 )
      {
         StX = X + nu * INCX;
 
         do
         {
            x0 = (*X);     x4 = X[incX4]; x1 = X[INCX ]; x5 = X[incX5];
            x2 = X[incX2]; x6 = X[incX6]; x3 = X[incX3]; x7 = X[incX7];
 
            *Y       = x0; Y[incY4] = x4; Y[INCY ] = x1; Y[incY5] = x5;
            Y[incY2] = x2; Y[incY6] = x6; Y[incY3] = x3; Y[incY7] = x7;
 
            X  += incX8;
            Y  += incY8;
 
         } while( X != StX );
      }
 
      for( i = N - nu; i != 0; i-- )
      {
         x0  = (*X);
         *Y  = x0;
 
         X  += INCX;
         Y  += INCY;
      }
   }
#endif
#ifdef HPL_CALL_FBLAS
#ifdef HPL_USE_F77_INTEGER_DEF
   const F77_INTEGER         F77N = N, F77incx = INCX, F77incy = INCY;
#else
#define F77N                 N
#define F77incx              INCX
#define F77incy              INCY
#endif
   F77dcopy( &F77N, X, &F77incx, Y, &F77incy );
#endif
/*
 * End of HPLAI_acopy
 */
}
 
#ifdef __cplusplus
}
#endif
