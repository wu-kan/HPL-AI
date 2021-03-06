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
#ifndef HPLAI_PANEL_HH
#define HPLAI_PANEL_HH
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hplai_misc.hh"

#include "hplai_pfact.hh"

#ifdef __cplusplus
extern "C"
{
#endif

    /*
 * ---------------------------------------------------------------------
 * panel function prototypes
 * ---------------------------------------------------------------------
 */

    void HPLAI_papanel_new
        STDC_ARGS((
            HPL_T_grid *,
            HPLAI_T_palg *,
            const int,
            const int,
            const int,
            HPLAI_T_pmat *,
            const int,
            const int,
            const int,
            HPLAI_T_panel **));
    void HPLAI_papanel_init
        STDC_ARGS((
            HPL_T_grid *,
            HPLAI_T_palg *,
            const int,
            const int,
            const int,
            HPLAI_T_pmat *,
            const int,
            const int,
            const int,
            HPLAI_T_panel *));
    int HPLAI_papanel_disp
        STDC_ARGS((
            HPLAI_T_panel **));
    int HPLAI_papanel_free
        STDC_ARGS((
            HPLAI_T_panel *));

#ifdef __cplusplus
}
#endif

#endif
/*
 * End of hplai_panel.hh
 */
