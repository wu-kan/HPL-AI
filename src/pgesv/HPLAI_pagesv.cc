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
    void HPLAI_pagesv(
        HPL_T_grid *GRID,
        HPLAI_T_palg *ALGO,
        HPLAI_T_pmat *A)
#else
void HPLAI_pagesv(GRID, ALGO, A)
    HPL_T_grid *GRID;
HPLAI_T_palg *ALGO;
HPLAI_T_pmat *A;
#endif
    {
        /* ..
 * .. Executable Statements ..
 */
        if (A->n <= 0)
            return;

        A->info = 0;

        if ((ALGO->depth == 0) || (GRID->npcol == 1))
        {
            HPLAI_pagesv0(GRID, ALGO, A);
        }
        else
        {
            HPLAI_pagesvK2(GRID, ALGO, A);
        }
        /*
 * Solve upper triangular system
 */
        if (A->info == 0)
            HPLAI_patrsv(GRID, A);
        /*
 * End of HPLAI_pagesv
 */
    }

#ifdef __cplusplus
}
#endif
