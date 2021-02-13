/*
 * Include files
 */
#include "hplai.h"

template <typename T>
static void HPLAI_aaxpy_template(
    const int N,
    const T ALPHA,
    const T *X,
    const int INCX,
    T *Y,
    const int INCY);

template <>
void HPLAI_aaxpy_template<double>(
    const int N,
    const double ALPHA,
    const double *X,
    const int INCX,
    double *Y,
    const int INCY)
{
   HPL_daxpy(N, ALPHA, X, INCX, Y, INCY);
}

#ifdef HPL_CALL_CBLAS
template <>
void HPLAI_aaxpy_template<float>(
    const int N,
    const float ALPHA,
    const float *X,
    const int INCX,
    float *Y,
    const int INCY)
{
   cblas_saxpy(N, ALPHA, X, INCX, Y, INCY);
}
#endif

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
      HPLAI_aaxpy_template(N, ALPHA, X, INCX, Y, INCY);
      /*
 * End of HPLAI_aaxpy
 */
   }

#ifdef __cplusplus
}
#endif
