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
#ifndef HPLAI_AUXIL_HH
#define HPLAI_AUXIL_HH
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#ifdef __cplusplus
extern "C"
{
#endif

#include "hpl_auxil.h"

#ifdef __cplusplus
}
#endif

#include "hplai_misc.hh"

#ifdef __cplusplus
extern "C"
{
#endif

#define HPLAI_NORM_A HPL_NORM_A
#define HPLAI_NORM_1 HPL_NORM_1
#define HPLAI_NORM_I HPL_NORM_I
#define HPLAI_T_NORM HPL_T_NORM

#define HPLAI_MACH_EPS HPL_MACH_EPS
#define HPLAI_MACH_SFMIN HPL_MACH_SFMIN
#define HPLAI_MACH_BASE HPL_MACH_BASE
#define HPLAI_MACH_PREC HPL_MACH_PREC
#define HPLAI_MACH_MLEN HPL_MACH_MLEN
#define HPLAI_MACH_RND HPL_MACH_RND
#define HPLAI_MACH_EMIN HPL_MACH_EMIN
#define HPLAI_MACH_RMIN HPL_MACH_RMIN
#define HPLAI_MACH_EMAX HPL_MACH_EMAX
#define HPLAI_MACH_RMAX HPL_MACH_RMAX
#define HPLAI_T_MACH HPL_T_MACH

#define HPLAI_fprintf HPL_fprintf
#define HPLAI_warn HPL_warn
#define HPLAI_abort HPL_abort

void HPLAI_alacpy
    STDC_ARGS((
        const int,
        const int,
        const HPLAI_T_AFLOAT *,
        const int,
        HPLAI_T_AFLOAT *,
        const int));
void HPLAI_alatcpy
    STDC_ARGS((
        const int,
        const int,
        const HPLAI_T_AFLOAT *,
        const int,
        HPLAI_T_AFLOAT *,
        const int));
HPLAI_T_AFLOAT HPLAI_alamch
    STDC_ARGS((
        const HPL_T_MACH));

#ifdef __cplusplus
}
#endif

#endif
/*
 * End of hplai_auxil.hh
 */
