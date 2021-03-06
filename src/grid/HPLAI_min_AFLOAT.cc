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

#ifdef STDC_HEADERS
    void HPLAI_min_AFLOAT(
        const int N,
        const void *IN,
        void *INOUT)
#else
void HPLAI_min_AFLOAT(N, IN, INOUT, DTYPE)
    const int N;
const void *IN;
void *INOUT
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_min_AFLOAT combines (min) two buffers.
 * 
 *
 * Arguments
 * =========
 *
 * N       (input)                       const int
 *         On entry, N  specifies  the  length  of  the  buffers  to  be
 *         combined. N must be at least zero.
 *
 * IN      (input)                       const void *
 *         On entry, IN points to the input-only buffer to be combined.
 *
 * INOUT   (input/output)                void *
 *         On entry, INOUT  points  to  the  input-output  buffer  to be
 *         combined.  On exit,  the  entries of this array contains  the
 *         combined results.
 *
 * DTYPE   (input)                       const HPL_T_TYPE
 *         On entry,  DTYPE  specifies the type of the buffers operands.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        register int i;
        {
            const HPLAI_T_AFLOAT *a = (const HPLAI_T_AFLOAT *)(IN);
            HPLAI_T_AFLOAT *b = (HPLAI_T_AFLOAT *)(INOUT);
            for (i = 0; i < N; i++)
                b[i] = Mmin(a[i], b[i]);
        }
        /*
 * End of HPLAI_min_AFLOAT
 */
    }

#ifdef __cplusplus
}
#endif
