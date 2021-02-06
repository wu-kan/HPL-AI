/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPLAI_aswap(
    const int N,
    HPLAI_T_AFLOAT *X,
    const int INCX,
    HPLAI_T_AFLOAT *Y,
    const int INCY)
#else
void HPLAI_aswap(N, X, INCX, Y, INCY)
    const int N;
HPLAI_T_AFLOAT *X;
const int INCX;
HPLAI_T_AFLOAT *Y;
const int INCY;
#endif
{
#ifdef HPL_CALL_CBLAS
   cblas_sswap(N, X, INCX, Y, INCY);
#endif
   /*
 * End of HPLAI_aswap
 */
}

#ifdef __cplusplus
}
#endif