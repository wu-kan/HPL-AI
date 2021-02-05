/*
 * Include files
 */
#include "hpl_ai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPL_AI_aswap(
    const int N,
    HPL_AI_T_afloat *X,
    const int INCX,
    HPL_AI_T_afloat *Y,
    const int INCY)
#else
void HPL_AI_aswap(N, X, INCX, Y, INCY)
    const int N;
HPL_AI_T_afloat *X;
const int INCX;
HPL_AI_T_afloat *Y;
const int INCY;
#endif
{
#ifdef HPL_CALL_CBLAS
   cblas_sswap(N, X, INCX, Y, INCY);
#endif
   /*
 * End of HPL_AI_aswap
 */
}

#ifdef __cplusplus
}
#endif