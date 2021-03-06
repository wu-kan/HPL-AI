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
#ifndef HPLAI_PTEST_H
#define HPLAI_PTEST_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_ptest.h"
/*
 * ---------------------------------------------------------------------
 * Data Structures
 * ---------------------------------------------------------------------
 */
#define HPLAI_S_test HPL_S_test
#define HPLAI_T_test HPL_T_test

/*
 * ---------------------------------------------------------------------
 * #define macro constants for testing only
 * ---------------------------------------------------------------------
 */
#define HPLAI_LINE_MAX HPL_LINE_MAX
#define HPLAI_MAX_PARAM HPL_MAX_PARAM
#define HPLAI_ISEED HPL_ISEED
/*
 * ---------------------------------------------------------------------
 * global timers for timing analysis only
 * ---------------------------------------------------------------------
 */
#ifdef HPL_DETAILED_TIMING
#define HPLAI_DETAILED_TIMING
#define HPLAI_TIMING_BEG HPL_TIMING_BEG
#define HPLAI_TIMING_N HPL_TIMING_N
#define HPLAI_TIMING_RPFACT HPL_TIMING_RPFACT
#define HPLAI_TIMING_PFACT HPL_TIMING_PFACT
#define HPLAI_TIMING_MXSWP HPL_TIMING_MXSWP
#define HPLAI_TIMING_UPDATE HPL_TIMING_UPDATE
#define HPLAI_TIMING_LASWP HPL_TIMING_LASWP
#define HPLAI_TIMING_PTRSV HPL_TIMING_PTRSV
#endif
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void HPLAI_pdinfo
    STDC_ARGS((
        HPLAI_T_test *,
        int *,
        int *,
        int *,
        int *,
        HPLAI_T_ORDER *,
        int *,
        int *,
        int *,
        int *,
        HPLAI_T_FACT *,
        int *,
        int *,
        int *,
        int *,
        int *,
        HPLAI_T_FACT *,
        int *,
        HPLAI_T_TOP *,
        int *,
        int *,
        HPLAI_T_SWAP *,
        int *,
        int *,
        int *,
        int *,
        int *));
void HPLAI_pdtest
    STDC_ARGS((
        HPLAI_T_test *,
        HPLAI_T_grid *,
        HPLAI_T_palg *,
        const int,
        const int));

#endif
/*
 * End of hplai_ptest.h
 */
