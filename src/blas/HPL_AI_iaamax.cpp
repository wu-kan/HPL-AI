/*
 * Include files
 */
#include "hpl_ai.h"

template <typename T>
static int HPL_AI_iaamax_template(
	const int N,
	const T *X,
	const int INCX);

template <>
int HPL_AI_iaamax_template<double>(
	const int N,
	const double *X,
	const int INCX)
{
	return HPL_idamax(N, X, INCX);
}

#ifdef HPL_CALL_CBLAS
template <>
int HPL_AI_iaamax_template<float>(
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
	int HPL_AI_iaamax(
		const int N,
		const HPL_AI_T_afloat *X,
		const int INCX)
#else
int HPL_idamax(N, X, INCX)
	const int N;
const HPL_AI_T_afloat *X;
const int INCX;
#endif
	{
		return HPL_AI_iaamax_template(N, X, INCX);
		/*
 * End of HPL_AI_iaamax
 */
	}

#ifdef __cplusplus
}
#endif
