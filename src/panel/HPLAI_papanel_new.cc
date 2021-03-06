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
    void HPLAI_papanel_new(
        HPLAI_T_grid *GRID,
        HPLAI_T_palg *ALGO,
        const int M,
        const int N,
        const int JB,
        HPLAI_T_pmat *A,
        const int IA,
        const int JA,
        const int TAG,
        HPLAI_T_panel **PANEL)
#else
void HPLAI_papanel_new(GRID, ALGO, M, N, JB, A, IA, JA, TAG, PANEL)
    HPLAI_T_grid *GRID;
HPL_T_palg *ALGO;
const int M;
const int N;
const int JB;
HPLAI_T_pmat *A;
const int IA;
const int JA;
const int TAG;
HPLAI_T_panel **PANEL;
#endif
    {
        /*
 * .. Local Variables ..
 */
        HPLAI_T_panel *p = NULL;
        /* ..
 * .. Executable Statements ..
 */
        /*
 * Allocate the panel structure - Check for enough memory
 */
        if (!(p = (HPLAI_T_panel *)malloc(sizeof(HPLAI_T_panel))))
        {
            HPL_pabort(__LINE__, "HPLAI_papanel_new", "Memory allocation failed");
        }

        HPLAI_papanel_init(GRID, ALGO, M, N, JB, A, IA, JA, TAG, p);
        *PANEL = p;
        /*
 * End of HPLAI_papanel_new
 */
    }

#ifdef __cplusplus
}
#endif
