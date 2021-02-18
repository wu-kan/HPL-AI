/*
 * Include files
 */
#include "hplai.h"

template <>
int64_t blas::iamax<double>(
	int64_t n,
	double const *x,
	int64_t incx)
{
	return HPL_idamax(n, x, incx);
}

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
	int HPLAI_iaamax(
		const int N,
		const HPLAI_T_AFLOAT *X,
		const int INCX)
#else
int HPL_idamax(N, X, INCX)
	const int N;
const HPLAI_T_AFLOAT *X;
const int INCX;
#endif
	{
		return blas::iamax(N, X, INCX);
		/*
 * End of HPLAI_iaamax
 */
	}

#ifdef __cplusplus
}
#endif
