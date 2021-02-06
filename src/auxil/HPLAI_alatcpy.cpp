/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
   void HPLAI_alatcpy(
       const int M,
       const int N,
       const HPLAI_T_AFLOAT *A,
       const int LDA,
       HPLAI_T_AFLOAT *B,
       const int LDB)
#else
void HPLAI_alatcpy(M, N, A, LDA, B, LDB)
    const int M;
const int N;
const HPLAI_T_AFLOAT *A;
const int LDA;
HPLAI_T_AFLOAT *B;
const int LDB;
#endif
   {
      /*
 * .. Local Variables ..
 */
      register int j;

      const HPLAI_T_AFLOAT *A0 = A;
      HPLAI_T_AFLOAT *B0 = B;

      /* ..
 * .. Executable Statements ..
 */
      if ((M <= 0) || (N <= 0))
         return;

      for (j = 0; j < N; j++, B0 += LDB)
         HPLAI_acopy(M, A0 + j, LDA, B0, 1);
      /*
 * End of HPLAI_alatcpy
 */
   }

#ifdef __cplusplus
}
#endif
