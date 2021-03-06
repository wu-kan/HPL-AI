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
#ifndef HPLAI_PAUXIL_HH
#define HPLAI_PAUXIL_HH
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hplai_misc.hh"
#include "hplai_auxil.hh"



#ifdef __cplusplus
extern "C"
{
#endif

#define HPLAI_indxg2lp HPL_indxg2lp
#define HPLAI_indxg2l HPL_indxg2l
#define HPLAI_indxg2p HPL_indxg2p
#define HPLAI_indxl2g HPL_indxl2g
#define HPLAI_infog2l HPL_infog2l
#define HPLAI_numroc HPL_numroc
#define HPLAI_numrocI HPL_numrocI

    void HPLAI_alaswp00N
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const int *));
    void HPLAI_alaswp10N
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const int *));
    void HPLAI_alaswp01N
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const int *,
            const int *));
    void HPLAI_alaswp01T
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const int *,
            const int *));
    void HPLAI_alaswp02N
        STDC_ARGS((
            const int,
            const int,
            const HPLAI_T_AFLOAT *,
            const int,
            HPLAI_T_AFLOAT *,
            HPLAI_T_AFLOAT *,
            const int,
            const int *,
            const int *));
    void HPLAI_alaswp03N
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const HPLAI_T_AFLOAT *,
            const HPLAI_T_AFLOAT *,
            const int));
    void HPLAI_alaswp03T
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const HPLAI_T_AFLOAT *,
            const HPLAI_T_AFLOAT *,
            const int));
    void HPLAI_alaswp04N
        STDC_ARGS((
            const int,
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const HPLAI_T_AFLOAT *,
            const HPLAI_T_AFLOAT *,
            const int,
            const int *,
            const int *));
    void HPLAI_alaswp04T
        STDC_ARGS((
            const int,
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const HPLAI_T_AFLOAT *,
            const HPLAI_T_AFLOAT *,
            const int,
            const int *,
            const int *));
    void HPLAI_alaswp05N
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const HPLAI_T_AFLOAT *,
            const int,
            const int *,
            const int *));
    void HPLAI_alaswp05T
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const HPLAI_T_AFLOAT *,
            const int,
            const int *,
            const int *));
    void HPLAI_alaswp06N
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const int *));
    void HPLAI_alaswp06T
        STDC_ARGS((
            const int,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            HPLAI_T_AFLOAT *,
            const int,
            const int *));

#define HPLAI_pabort HPL_pabort
#define HPLAI_pwarn HPL_pwarn

#ifdef __cplusplus
}
#endif

#endif
/*
 * End of hplai_pauxil.hh
 */
