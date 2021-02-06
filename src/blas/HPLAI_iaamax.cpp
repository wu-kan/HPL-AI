/*
 * Include files
 */
#include "hplai.h"

template <typename T>
static int HPLAI_iaamax_template(
	const int N,
	const T *X,
	const int INCX);

template <>
int HPLAI_iaamax_template<double>(
	const int N,
	const double *X,
	const int INCX)
{
	return HPL_idamax(N, X, INCX);
}

#ifdef HPL_CALL_CBLAS
template <>
int HPLAI_iaamax_template<float>(
	const int N,
	const float *X,
	const int INCX)
{
	return ((int)(cblas_isamax(N, X, INCX)));
}
#endif

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
		return HPLAI_iaamax_template(N, X, INCX);
		/*
 * End of HPLAI_iaamax
 */
	}

#ifdef __cplusplus
}
#endif
