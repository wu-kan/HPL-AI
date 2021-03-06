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
    void HPLAI_copyL(
        HPLAI_T_panel *PANEL)
#else
void HPLAI_copyL(PANEL)
    HPLAI_T_panel *PANEL;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_copyL copies  the  panel of columns, the L1 replicated submatrix,
 * the pivot array  and  the info scalar into a contiguous workspace for
 * later broadcast.
 *  
 * The copy of this panel  into  a contiguous buffer  can be enforced by
 * specifying -DHPL_COPY_L in the architecture specific Makefile.
 *
 * Arguments
 * =========
 *
 * PANEL   (input/output)                HPLAI_T_panel *
 *         On entry,  PANEL  points to the  current panel data structure
 *         being broadcast.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        int jb, lda;
        /* ..
 * .. Executable Statements ..
 */
        if (PANEL->grid->mycol == PANEL->pcol)
        {
            jb = PANEL->jb;
            lda = PANEL->lda;

            if (PANEL->grid->myrow == PANEL->prow)
            {
                HPLAI_alacpy(PANEL->mp - jb, jb, Mptr(PANEL->A, jb, -jb, lda),
                             lda, PANEL->L2, PANEL->ldl2);
            }
            else
            {
                HPLAI_alacpy(PANEL->mp, jb, Mptr(PANEL->A, 0, -jb, lda),
                             lda, PANEL->L2, PANEL->ldl2);
            }
        }
        /*
 * End of HPLAI_copyL
 */
    }

#ifdef __cplusplus
}
#endif
