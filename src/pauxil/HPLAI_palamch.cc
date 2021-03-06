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
    HPLAI_T_AFLOAT HPLAI_palamch(
        MPI_Comm COMM,
        const HPLAI_T_MACH CMACH)
#else
HPLAI_T_AFLOAT HPLAI_palamch(COMM, CMACH)
    MPI_Comm COMM;
const HPLAI_T_MACH CMACH;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_palamch determines  machine-specific  arithmetic  constants  such  as
 * the relative machine precision (eps),  the safe minimum(sfmin) such that
 * 1/sfmin does not overflow, the base of the machine (base), the precision
 * (prec),  the  number  of  (base)  digits in the  mantissa  (t),  whether
 * rounding occurs in addition (rnd = 1.0 and 0.0 otherwise),  the  minimum
 * exponent before  (gradual)  underflow (emin),  the  underflow  threshold
 * (rmin)- base**(emin-1), the largest exponent before overflow (emax), the
 * overflow threshold (rmax)  - (base**emax)*(1-eps).
 *
 * Arguments
 * =========
 *
 * COMM    (global/local input)          MPI_Comm
 *         The MPI communicator identifying the process collection.
 *
 * CMACH   (global input)                const HPLAI_T_MACH
 *         Specifies the value to be returned by HPLAI_palamch            
 *            = HPLAI_MACH_EPS,   HPLAI_palamch := eps (default)            
 *            = HPLAI_MACH_SFMIN, HPLAI_palamch := sfmin                    
 *            = HPLAI_MACH_BASE,  HPLAI_palamch := base                     
 *            = HPLAI_MACH_PREC,  HPLAI_palamch := eps*base                 
 *            = HPLAI_MACH_MLEN,  HPLAI_palamch := t                        
 *            = HPLAI_MACH_RND,   HPLAI_palamch := rnd                      
 *            = HPLAI_MACH_EMIN,  HPLAI_palamch := emin                     
 *            = HPLAI_MACH_RMIN,  HPLAI_palamch := rmin                     
 *            = HPLAI_MACH_EMAX,  HPLAI_palamch := emax                     
 *            = HPLAI_MACH_RMAX,  HPLAI_palamch := rmax                     
 *          
 *         where                                                        
 *          
 *            eps   = relative machine precision,                       
 *            sfmin = safe minimum,                                     
 *            base  = base of the machine,                              
 *            prec  = eps*base,                                         
 *            t     = number of digits in the mantissa,                 
 *            rnd   = 1.0 if rounding occurs in addition,               
 *            emin  = minimum exponent before underflow,                
 *            rmin  = underflow threshold,                              
 *            emax  = largest exponent before overflow,                 
 *            rmax  = overflow threshold.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        HPLAI_T_AFLOAT param;
        /* ..
 * .. Executable Statements ..
 */
        param = HPLAI_alamch(CMACH);

        switch (CMACH)
        {
        case HPLAI_MACH_EPS:
        case HPLAI_MACH_SFMIN:
        case HPLAI_MACH_EMIN:
        case HPLAI_MACH_RMIN:
            (void)HPLAI_all_reduce_AFLOAT((void *)(&param), 1,
                                          HPLAI_max_AFLOAT, COMM);
            break;
        case HPLAI_MACH_EMAX:
        case HPLAI_MACH_RMAX:
            (void)HPLAI_all_reduce_AFLOAT((void *)(&param), 1,
                                          HPLAI_min_AFLOAT, COMM);
            break;
        default:
            break;
        }

        return (param);
        /*
 * End of HPLAI_palamch
 */
    }

#ifdef __cplusplus
}
#endif
