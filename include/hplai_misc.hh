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
#ifndef HPLAI_MISC_HH
#define HPLAI_MISC_HH
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#ifdef __cplusplus
extern "C"
{
#endif

#include "hpl.h"

#ifdef __cplusplus
}
#endif

#include <math.h>
#include <float.h>
//use blaspp https://bitbucket.org/icl/blaspp/src/master/
#include <blas.hh>

/*
 * ---------------------------------------------------------------------
 * #define macros definitions
 * ---------------------------------------------------------------------
 */
#ifndef HPLAI_T_AFLOAT
#define HPLAI_T_AFLOAT float
#endif

#define HPLAI_rone ((HPLAI_T_AFLOAT)HPL_rone)
#define HPLAI_rtwo ((HPLAI_T_AFLOAT)HPL_rtwo)
#define HPLAI_rzero ((HPLAI_T_AFLOAT)HPL_rzero)

#define HPLAI_PTR HPL_PTR

#endif
/*
 * End of hplai_misc.hh
 */
