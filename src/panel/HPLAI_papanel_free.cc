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
#include "hplai.hh"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    int HPLAI_papanel_free(
        HPLAI_T_panel *PANEL)
#else
int HPLAI_papanel_free(PANEL)
    HPLAI_T_panel *PANEL;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPL_pdpanel_free deallocates  the panel resources  and  stores the error
 * code returned by the panel factorization.
 *
 * Arguments
 * =========
 *
 * PANEL   (local input/output)          HPL_T_panel *
 *         On entry,  PANEL  points  to  the  panel data  structure from
 *         which the resources should be deallocated.
 *
 * ---------------------------------------------------------------------
 */
        /* ..
 * .. Executable Statements ..
 */
        if (PANEL->pmat->info == 0)
            PANEL->pmat->info = *(PANEL->DINFO);
#ifdef HPL_CALL_VSIPL
        /*
 * Release the blocks
 */
        (void)vsip_blockrelease_d(PANEL->L1block, VSIP_TRUE);
        (void)vsip_blockrelease_d(PANEL->L2block, VSIP_TRUE);
        if (PANEL->grid->nprow > 1)
            (void)vsip_blockrelease_d(PANEL->Ublock, VSIP_TRUE);
        /*
 * Destroy blocks
 */
        vsip_blockdestroy_d(PANEL->L1block);
        vsip_blockdestroy_d(PANEL->L2block);
        if (PANEL->grid->nprow > 1)
            vsip_blockdestroy_d(PANEL->Ublock);
#endif

        if (PANEL->WORK)
            free(PANEL->WORK);
        if (PANEL->IWORK)
            free(PANEL->IWORK);

        return (MPI_SUCCESS);
        /*
 * End of HPLAI_papanel_free
 */
    }

#ifdef __cplusplus
}
#endif
