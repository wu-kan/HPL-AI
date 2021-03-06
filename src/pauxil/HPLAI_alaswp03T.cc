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
#ifndef HPL_LASWP03T_DEPTH
#define HPL_LASWP03T_DEPTH 32
#define HPL_LASWP03T_LOG2_DEPTH 5
#endif

#ifdef STDC_HEADERS
    void HPLAI_alaswp03T(
        const int M,
        const int N,
        HPLAI_T_AFLOAT *U,
        const int LDU,
        const HPLAI_T_AFLOAT *W0,
        const HPLAI_T_AFLOAT *W,
        const int LDW)
#else
void HPLAI_alaswp03T(M, N, U, LDU, W0, W, LDW)
    const int M;
const int N;
HPLAI_T_AFLOAT *U;
const int LDU;
const HPLAI_T_AFLOAT *W0;
const HPLAI_T_AFLOAT *W;
const int LDW;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_alaswp03T copies  columns of W into an array U.  The  destination
 * in U of these columns contained in W is stored within W0.
 *
 * Arguments
 * =========
 *
 * M       (local input)                 const int
 *         On entry, M  specifies  the  number  of columns of  W  stored
 *         contiguously that should be copied into U. M must be at least
 *         zero.
 *
 * N       (local input)                 const int
 *         On entry,  N  specifies  the  length of columns of  W  stored
 *         contiguously that should be copied into U. N must be at least
 *         zero.
 *
 * U       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry, U points to an array of dimension (LDU,M).  Columns
 *         of W are copied within the array U at the positions specified
 *         in W0.
 *
 * LDU     (local input)                 const int
 *         On entry, LDU specifies the leading dimension of the array U.
 *         LDU must be at least MAX(1,N).
 *
 * W0      (local input)                 const HPLAI_T_AFLOAT *
 *         On entry,  W0  is an array of size (M-1)*LDW+1, that contains
 *         the destination offset  in U where the columns of W should be
 *         copied.
 *
 * W       (local input)                 const HPLAI_T_AFLOAT *
 *         On entry, W  is an array of size (LDW,M),  that contains data
 *         to be copied into U. For i in [0..M),  entries W(:,i)  should
 *         be copied into the row or column W0(i*LDW) of U.
 *
 * LDW     (local input)                 const int
 *         On entry, LDW specifies the leading dimension of the array W.
 *         LDW must be at least MAX(1,N+1).
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        const HPLAI_T_AFLOAT *w = W, *w0;
        HPLAI_T_AFLOAT *u0;
        const int incU = (1 << HPL_LASWP03T_LOG2_DEPTH);
        int nr, nu;
        register int i, j;
        /* ..
 * .. Executable Statements ..
 */
        if ((M <= 0) || (N <= 0))
            return;

        nr = N - (nu = (int)(((unsigned int)(N) >> HPL_LASWP03T_LOG2_DEPTH) << HPL_LASWP03T_LOG2_DEPTH));

        for (j = 0; j < nu;
             j += HPL_LASWP03T_DEPTH, U += incU, w += HPL_LASWP03T_DEPTH)
        {
            for (i = 0; i < M; i++)
            {
                u0 = U + (size_t)(*(W0 + (size_t)(i) * (size_t)(LDW))) * (size_t)(LDU);
                w0 = w + (size_t)(i) * (size_t)(LDW);

                u0[0] = w0[0];
#if (HPL_LASWP03T_DEPTH > 1)
                u0[1] = w0[1];
#endif
#if (HPL_LASWP03T_DEPTH > 2)
                u0[2] = w0[2];
                u0[3] = w0[3];
#endif
#if (HPL_LASWP03T_DEPTH > 4)
                u0[4] = w0[4];
                u0[5] = w0[5];
                u0[6] = w0[6];
                u0[7] = w0[7];
#endif
#if (HPL_LASWP03T_DEPTH > 8)
                u0[8] = w0[8];
                u0[9] = w0[9];
                u0[10] = w0[10];
                u0[11] = w0[11];
                u0[12] = w0[12];
                u0[13] = w0[13];
                u0[14] = w0[14];
                u0[15] = w0[15];
#endif
#if (HPL_LASWP03T_DEPTH > 16)
                u0[16] = w0[16];
                u0[17] = w0[17];
                u0[18] = w0[18];
                u0[19] = w0[19];
                u0[20] = w0[20];
                u0[21] = w0[21];
                u0[22] = w0[22];
                u0[23] = w0[23];
                u0[24] = w0[24];
                u0[25] = w0[25];
                u0[26] = w0[26];
                u0[27] = w0[27];
                u0[28] = w0[28];
                u0[29] = w0[29];
                u0[30] = w0[30];
                u0[31] = w0[31];
#endif
            }
        }

        if (nr > 0)
        {
            for (i = 0; i < M; i++)
            {
                u0 = U + (size_t)(*(W0 + (size_t)(i) * (size_t)(LDW))) * (size_t)(LDU);
                w0 = w + (size_t)(i) * (size_t)(LDW);
                for (j = 0; j < nr; j++)
                {
                    u0[j] = w0[j];
                }
            }
        }
        /*
 * End of HPLAI_alaswp03T
 */
    }

#ifdef __cplusplus
}
#endif
