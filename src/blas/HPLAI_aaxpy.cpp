/*
 * Include files
 */
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

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
   void HPLAI_aaxpy(
       const int N,
       const HPLAI_T_AFLOAT ALPHA,
       const HPLAI_T_AFLOAT *X,
       const int INCX,
       HPLAI_T_AFLOAT *Y,
       const int INCY)
#else
void HPLAI_aaxpy(N, ALPHA, X, INCX, Y, INCY)
    const int N;
const HPLAI_T_AFLOAT ALPHA;
const HPLAI_T_AFLOAT *X;
const int INCX;
HPLAI_T_AFLOAT *Y;
const int INCY;
#endif
   {
      blas::axpy(N, ALPHA, X, INCX, Y, INCY);
      /*
 * End of HPLAI_aaxpy
 */
   }

#ifdef __cplusplus
}
#endif
