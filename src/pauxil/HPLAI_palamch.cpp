/* 
 * -- High Performance Computing Linpack Benchmark (HPL)                
 *    HPL - 2.3 - December 2, 2018                          
 *    Antoine P. Petitet                                                
 *    University of Tennessee, Knoxville                                
 *    Innovative Computing Laboratory                                 
 *    (C) Copyright 2000-2008 All Rights Reserved                       
 *                                                                      
 * -- Copyright notice and Licensing terms:                             
 *                                                                      
 * Redistribution  and  use in  source and binary forms, with or without
 * modification, are  permitted provided  that the following  conditions
 * are met:                                                             
 *                                                                      
 * 1. Redistributions  of  source  code  must retain the above copyright
 * notice, this list of conditions and the following disclaimer.        
 *                                                                      
 * 2. Redistributions in binary form must reproduce  the above copyright
 * notice, this list of conditions,  and the following disclaimer in the
 * documentation and/or other materials provided with the distribution. 
 *                                                                      
 * 3. All  advertising  materials  mentioning  features  or  use of this
 * software must display the following acknowledgement:                 
 * This  product  includes  software  developed  at  the  University  of
 * Tennessee, Knoxville, Innovative Computing Laboratory.             
 *                                                                      
 * 4. The name of the  University,  the name of the  Laboratory,  or the
 * names  of  its  contributors  may  not  be used to endorse or promote
 * products  derived   from   this  software  without  specific  written
 * permission.                                                          
 *                                                                      
 * -- Disclaimer:                                                       
 *                                                                      
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR  CONTRIBUTORS  BE  LIABLE FOR ANY  DIRECT,  INDIRECT,  INCIDENTAL,
 * SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES  (INCLUDING,  BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,  OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 * ---------------------------------------------------------------------
 */ 
/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
HPLAI_T_AFLOAT HPLAI_palamch
(
   MPI_Comm                         COMM,
   const HPLAI_T_MACH                 CMACH
)
#else
HPLAI_T_AFLOAT HPLAI_palamch
( COMM, CMACH )
   MPI_Comm                         COMM;
   const HPLAI_T_MACH                 CMACH;
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
   HPLAI_T_AFLOAT                     param;
/* ..
 * .. Executable Statements ..
 */
   param = HPLAI_alamch( CMACH );

   switch( CMACH )
   {
      case HPLAI_MACH_EPS   :
      case HPLAI_MACH_SFMIN :
      case HPLAI_MACH_EMIN  :
      case HPLAI_MACH_RMIN  :
         (void) HPLAI_all_reduce_AFLOAT( (void *)(&param), 1,
                                HPLAI_max_AFLOAT, COMM );
         break;
      case HPLAI_MACH_EMAX  :
      case HPLAI_MACH_RMAX  :
         (void) HPLAI_all_reduce_AFLOAT( (void *)(&param), 1,
                                HPLAI_min_AFLOAT, COMM );
         break;
      default             :
         break;
   } 

   return( param );
/*
 * End of HPLAI_palamch
 */
}

#ifdef __cplusplus
}
#endif
