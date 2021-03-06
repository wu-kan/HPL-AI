/*
 * MIT License
 * 
 * Copyright (c) 2021 WuK
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
/*
 * Include files
 */
#include "hplai.hh"

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * Define default value for unrolling factor
 */
#ifndef HPL_LASWP01N_DEPTH
#define HPL_LASWP01N_DEPTH 32
#define HPL_LASWP01N_LOG2_DEPTH 5
#endif

#ifdef STDC_HEADERS
    void HPLAI_alaswp01N(
        const int M,
        const int N,
        HPLAI_T_AFLOAT *A,
        const int LDA,
        HPLAI_T_AFLOAT *U,
        const int LDU,
        const int *LINDXA,
        const int *LINDXAU)
#else
void HPLAI_alaswp01N(M, N, A, LDA, U, LDU, LINDXA, LINDXAU)
    const int M;
const int N;
HPLAI_T_AFLOAT *A;
const int LDA;
HPLAI_T_AFLOAT *U;
const int LDU;
const int *LINDXA;
const int *LINDXAU;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_alaswp01N copies  scattered rows  of  A  into itself  and into an
 * array  U.  The row offsets in  A  of the source rows are specified by
 * LINDXA.  The  destination of those rows are specified by  LINDXAU.  A
 * positive value of  LINDXAU indicates that the array destination is U,
 * and A otherwise.
 *
 * Arguments
 * =========
 *
 * M       (local input)                 const int
 *         On entry, M  specifies the number of rows of A that should be
 *         moved within A or copied into U. M must be at least zero.
 *
 * N       (local input)                 const int
 *         On entry, N  specifies the length of rows of A that should be
 *         moved within A or copied into U. N must be at least zero.
 *
 * A       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry, A points to an array of dimension (LDA,N). The rows
 *         of this array specified by LINDXA should be moved within A or
 *         copied into U.
 *
 * LDA     (local input)                 const int
 *         On entry, LDA specifies the leading dimension of the array A.
 *         LDA must be at least MAX(1,M).
 *
 * U       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry, U points to an array of dimension (LDU,N). The rows
 *         of A specified by LINDXA are be copied within this array U at
 *         the positions indicated by positive values of LINDXAU.
 *
 * LDU     (local input)                 const int
 *         On entry, LDU specifies the leading dimension of the array U.
 *         LDU must be at least MAX(1,M).
 *
 * LINDXA  (local input)                 const int *
 *         On entry, LINDXA is an array of dimension M that contains the
 *         local  row indexes  of  A  that should be moved within  A  or
 *         or copied into U.
 *
 * LINDXAU (local input)                 const int *
 *         On entry, LINDXAU  is an array of dimension  M that  contains
 *         the local  row indexes of  U  where the rows of  A  should be
 *         copied at. This array also contains the  local row offsets in
 *         A where some of the rows of A should be moved to.  A positive
 *         value of  LINDXAU[i]  indicates that the row  LINDXA[i]  of A
 *         should be copied into U at the position LINDXAU[i]; otherwise
 *         the row  LINDXA[i]  of  A  should be moved  at  the  position
 *         -LINDXAU[i] within A.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        HPLAI_T_AFLOAT *a0, *a1;
        const int incA = (int)((unsigned int)(LDA) << HPL_LASWP01N_LOG2_DEPTH),
                  incU = (int)((unsigned int)(LDU) << HPL_LASWP01N_LOG2_DEPTH);
        int lda1, nu, nr;
        register int i, j;
        /* ..
 * .. Executable Statements ..
 */
        if ((M <= 0) || (N <= 0))
            return;

        nr = N - (nu = (int)(((unsigned int)(N) >> HPL_LASWP01N_LOG2_DEPTH) << HPL_LASWP01N_LOG2_DEPTH));

        for (j = 0; j < nu; j += HPL_LASWP01N_DEPTH, A += incA, U += incU)
        {
            for (i = 0; i < M; i++)
            {
                a0 = A + (size_t)(LINDXA[i]);
                if (LINDXAU[i] >= 0)
                {
                    a1 = U + (size_t)(LINDXAU[i]);
                    lda1 = LDU;
                }
                else
                {
                    a1 = A - (size_t)(LINDXAU[i]);
                    lda1 = LDA;
                }

                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
#if (HPL_LASWP01N_DEPTH > 1)
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
#endif
#if (HPL_LASWP01N_DEPTH > 2)
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
#endif
#if (HPL_LASWP01N_DEPTH > 4)
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
#endif
#if (HPL_LASWP01N_DEPTH > 8)
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
#endif
#if (HPL_LASWP01N_DEPTH > 16)
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
                *a1 = *a0;
                a1 += lda1;
                a0 += LDA;
#endif
            }
        }

        if (nr)
        {
            for (i = 0; i < M; i++)
            {
                a0 = A + (size_t)(LINDXA[i]);
                if (LINDXAU[i] >= 0)
                {
                    a1 = U + (size_t)(LINDXAU[i]);
                    lda1 = LDU;
                }
                else
                {
                    a1 = A - (size_t)(LINDXAU[i]);
                    lda1 = LDA;
                }
                for (j = 0; j < nr; j++, a1 += lda1, a0 += LDA)
                {
                    *a1 = *a0;
                }
            }
        }
        /*
 * End of HPLAI_alaswp01N
 */
    }

#ifdef __cplusplus
}
#endif
