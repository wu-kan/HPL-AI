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
#ifndef HPL_LASWP06N_DEPTH
#define HPL_LASWP06N_DEPTH 32
#define HPL_LASWP06N_LOG2_DEPTH 5
#endif

#ifdef STDC_HEADERS
    void HPLAI_alaswp06N(
        const int M,
        const int N,
        HPLAI_T_AFLOAT *A,
        const int LDA,
        HPLAI_T_AFLOAT *U,
        const int LDU,
        const int *LINDXA)
#else
void HPLAI_alaswp06N(M, N, A, LDA, U, LDU, LINDXA)
    const int M;
const int N;
HPLAI_T_AFLOAT *A;
const int LDA;
HPLAI_T_AFLOAT *U;
const int LDU;
const int *LINDXA;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_alaswp06N swaps rows of  U  with rows of A at positions
 * indicated by LINDXA.
 *
 * Arguments
 * =========
 *
 * M       (local input)                 const int
 *         On entry, M  specifies the number of rows of A that should be
 *         swapped with rows of U. M must be at least zero.
 *
 * N       (local input)                 const int
 *         On entry, N specifies the length of the rows of A that should
 *         be swapped with rows of U. N must be at least zero.
 *
 * A       (local output)                HPLAI_T_AFLOAT *
 *         On entry, A points to an array of dimension (LDA,N). On exit,
 *         the  rows of this array specified by  LINDXA  are replaced by
 *         rows or columns of U.
 *
 * LDA     (local input)                 const int
 *         On entry, LDA specifies the leading dimension of the array A.
 *         LDA must be at least MAX(1,M).
 *
 * U       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry,  U  points  to an array of dimension (LDU,N).  This
 *         array contains the rows of U that are to be swapped with rows
 *         of A.
 *
 * LDU     (local input)                 const int
 *         On entry, LDU specifies the leading dimension of the array U.
 *         LDU must be at least MAX(1,M).
 *
 * LINDXA  (local input)                 const int *
 *         On entry, LINDXA is an array of dimension M that contains the
 *         local row indexes of A that should be swapped with U.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        HPLAI_T_AFLOAT r;
        HPLAI_T_AFLOAT *U0 = U, *a0, *u0;
        const int incA = (int)((unsigned int)(LDA) << HPL_LASWP06N_LOG2_DEPTH),
                  incU = (int)((unsigned int)(LDU) << HPL_LASWP06N_LOG2_DEPTH);
        int nr, nu;
        register int i, j;
        /* ..
 * .. Executable Statements ..
 */
        if ((M <= 0) || (N <= 0))
            return;

        nr = N - (nu = (int)(((unsigned int)(N) >> HPL_LASWP06N_LOG2_DEPTH) << HPL_LASWP06N_LOG2_DEPTH));

        for (j = 0; j < nu; j += HPL_LASWP06N_DEPTH, A += incA, U0 += incU)
        {
            for (i = 0; i < M; i++)
            {
                a0 = A + (size_t)(LINDXA[i]);
                u0 = U0 + (size_t)(i);

                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
#if (HPL_LASWP06N_DEPTH > 1)
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
#endif
#if (HPL_LASWP06N_DEPTH > 2)
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
#endif
#if (HPL_LASWP06N_DEPTH > 4)
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
#endif
#if (HPL_LASWP06N_DEPTH > 8)
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
#endif
#if (HPL_LASWP06N_DEPTH > 16)
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
                r = *a0;
                *a0 = *u0;
                *u0 = r;
                a0 += LDA;
                u0 += LDU;
#endif
            }
        }

        if (nr)
        {
            for (i = 0; i < M; i++)
            {
                a0 = A + (size_t)(LINDXA[i]);
                u0 = U0 + (size_t)(i);
                for (j = 0; j < nr; j++, a0 += LDA, u0 += LDU)
                {
                    r = *a0;
                    *a0 = *u0;
                    *u0 = r;
                }
            }
        }
        /*
 * End of HPLAI_alaswp06N
 */
    }

#ifdef __cplusplus
}
#endif
