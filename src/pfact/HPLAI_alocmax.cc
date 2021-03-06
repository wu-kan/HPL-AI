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
    void HPLAI_alocmax(
        HPLAI_T_panel *PANEL,
        const int N,
        const int II,
        const int JJ,
        HPLAI_T_AFLOAT *WORK)
#else
void HPLAI_alocmax(PANEL, N, II, JJ, WORK)
    HPLAI_T_panel *PANEL;
const int N;
const int II;
const int JJ;
HPLAI_T_AFLOAT *WORK;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_alocmax finds  the maximum entry in the current column  and packs
 * the useful information in  WORK[0:3].  On exit,  WORK[0] contains the
 * local maximum  absolute value  scalar,  WORK[1] is the  corresponding
 * local row index,  WORK[2]  is the corresponding global row index, and
 * WORK[3] is the coordinate of the process owning this max.  When N  is
 * less than 1, the WORK[0:2] is initialized to zero, and WORK[3] is set
 * to the total number of process rows.
 *
 * Arguments
 * =========
 *
 * PANEL   (local input/output)          HPLAI_T_panel *
 *         On entry,  PANEL  points to the data structure containing the
 *         panel information.
 *
 * N       (local input)                 const int
 *         On entry,  N specifies the local number of rows of the column
 *         of A on which we operate.
 *
 * II      (local input)                 const int
 *         On entry, II  specifies the row offset where the column to be
 *         operated on starts with respect to the panel.
 *
 * JJ      (local input)                 const int
 *         On entry, JJ  specifies the column offset where the column to
 *         be operated on starts with respect to the panel.
 *
 * WORK    (local workspace)             HPLAI_T_AFLOAT *
 *         On entry, WORK  is  a workarray of size at least 4.  On exit,
 *         WORK[0] contains  the  local  maximum  absolute value scalar,
 *         WORK[1] contains  the corresponding local row index,  WORK[2]
 *         contains the corresponding global row index, and  WORK[3]  is
 *         the coordinate of process owning this max.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        HPLAI_T_AFLOAT *A;
        int kk, igindx, ilindx, myrow, nb, nprow;
        /* ..
 * .. Executable Statements ..
 */
        if (N > 0)
        {
            A = Mptr(PANEL->A, II, JJ, PANEL->lda);
            myrow = PANEL->grid->myrow;
            nprow = PANEL->grid->nprow;
            nb = PANEL->nb;
            kk = PANEL->ii + II + (ilindx = blas::iamax<HPLAI_T_AFLOAT>(N, A, 1));
            Mindxl2g(igindx, kk, nb, nb, myrow, 0, nprow);
            /*
 * WORK[0] := local maximum absolute value scalar,
 * WORK[1] := corresponding local  row index,
 * WORK[2] := corresponding global row index,
 * WORK[3] := coordinate of process owning this max.
 */
            WORK[0] = A[ilindx];
            WORK[1] = (HPLAI_T_AFLOAT)(ilindx);
            WORK[2] = (HPLAI_T_AFLOAT)(igindx);
            WORK[3] = (HPLAI_T_AFLOAT)(myrow);
        }
        else
        {
            /*
 * If I do not have any row of A, then set the coordinate of the process
 * (WORK[3]) owning this "ghost" row,  such that it  will never be used,
 * even if there are only zeros in the current column of A.
 */
            WORK[0] = WORK[1] = WORK[2] = HPLAI_rzero;
            WORK[3] = (HPLAI_T_AFLOAT)(PANEL->grid->nprow);
        }
        /*
 * End of HPLAI_alocmax
 */
    }

#ifdef __cplusplus
}
#endif
