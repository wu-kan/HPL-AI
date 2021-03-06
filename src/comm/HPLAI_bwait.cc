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
    int HPLAI_bwait(
        HPLAI_T_panel *PANEL)
#else
int HPLAI_bwait(PANEL)
    HPLAI_T_panel *PANEL;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_bwait HPLAI_bwait waits  for  the  row  broadcast  of  the current  panel  to
 * terminate.  Successful completion is indicated by the returned  error
 * code HPL_SUCCESS.
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
        int ierr;
        HPL_T_TOP top;
        /* ..
 * .. Executable Statements ..
 */
        if (PANEL->grid->npcol <= 1)
            return (HPL_SUCCESS);
        /*
 * Retrieve the selected virtual broadcast topology
 */
        top = PANEL->algo->btopo;

        switch (top)
        {
        case HPL_1RING_M:
            ierr = HPLAI_bwait_1rinM(PANEL);
            break;
        case HPL_1RING:
            ierr = HPLAI_bwait_1ring(PANEL);
            break;
        case HPL_2RING_M:
            ierr = HPLAI_bwait_2rinM(PANEL);
            break;
        case HPL_2RING:
            ierr = HPLAI_bwait_2ring(PANEL);
            break;
        case HPL_BLONG_M:
            ierr = HPLAI_bwait_blonM(PANEL);
            break;
        case HPL_BLONG:
            ierr = HPLAI_bwait_blong(PANEL);
            break;
        default:
            ierr = HPL_SUCCESS;
        }

        return (ierr);
        /*
 * End of HPLAI_bwait
 */
    }

#ifdef __cplusplus
}
#endif
