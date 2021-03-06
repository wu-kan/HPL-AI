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
    void HPLAI_pafact(
        HPLAI_T_panel *PANEL)
#else
void HPLAI_pafact(PANEL)
    HPLAI_T_panel *PANEL;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_pafact recursively factorizes a  1-dimensional  panel of columns.
 * The  RPFACT  function pointer specifies the recursive algorithm to be
 * used, either Crout, Left- or Right looking.  NBMIN allows to vary the
 * recursive stopping criterium in terms of the number of columns in the
 * panel, and  NDIV allows to specify the number of subpanels each panel
 * should be divided into. Usuallly a value of 2 will be chosen. Finally
 * PFACT is a function pointer specifying the non-recursive algorithm to
 * to be used on at most NBMIN columns. One can also choose here between
 * Crout, Left- or Right looking.  Empirical tests seem to indicate that
 * values of 4 or 8 for NBMIN give the best results.
 *  
 * Bi-directional  exchange  is  used  to  perform  the  swap::broadcast
 * operations  at once  for one column in the panel.  This  results in a
 * lower number of slightly larger  messages than usual.  On P processes
 * and assuming bi-directional links,  the running time of this function
 * can be approximated by (when N is equal to N0):                      
 *  
 *    N0 * log_2( P ) * ( lat + ( 2*N0 + 4 ) / bdwth ) +
 *    N0^2 * ( M - N0/3 ) * gam2-3
 *  
 * where M is the local number of rows of  the panel, lat and bdwth  are
 * the latency and bandwidth of the network for  HPLAI_T_AFLOAT  precision  real
 * words, and  gam2-3  is  an estimate of the  Level 2 and Level 3  BLAS
 * rate of execution. The  recursive  algorithm  allows indeed to almost
 * achieve  Level 3 BLAS  performance  in the panel factorization.  On a
 * large  number of modern machines,  this  operation is however latency
 * bound,  meaning  that its cost can  be estimated  by only the latency
 * portion N0 * log_2(P) * lat.  Mono-directional links will HPLAI_T_AFLOAT this
 * communication cost.
 *
 * Arguments
 * =========
 *
 * PANEL   (local input/output)          HPLAI_T_panel *
 *         On entry,  PANEL  points to the data structure containing the
 *         panel information.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        void *vptr = NULL;
        int align, jb;
        /* ..
 * .. Executable Statements ..
 */
        jb = PANEL->jb;
        PANEL->n -= jb;
        PANEL->ja += jb;

        if ((PANEL->grid->mycol != PANEL->pcol) || (jb <= 0))
            return;
#ifdef HPL_DETAILED_TIMING
        HPL_ptimer(HPL_TIMING_RPFACT);
#endif
        align = PANEL->algo->align;
        vptr = (void *)malloc(((size_t)(align) +
                               (size_t)(((4 + ((unsigned int)(jb) << 1)) << 1))) *
                              sizeof(HPLAI_T_AFLOAT));
        if (vptr == NULL)
        {
            HPLAI_pabort(__LINE__, "HPLAI_pafact", "Memory allocation failed");
        }
        /*
 * Factor the panel - Update the panel pointers
 */
        PANEL->algo->rffun(PANEL, PANEL->mp, jb, 0, (HPLAI_T_AFLOAT *)HPL_PTR(vptr, ((size_t)(align) * sizeof(HPLAI_T_AFLOAT))));
        if (vptr)
            free(vptr);

        PANEL->A = Mptr(PANEL->A, 0, jb, PANEL->lda);
        PANEL->nq -= jb;
        PANEL->jj += jb;
#ifdef HPL_DETAILED_TIMING
        HPL_ptimer(HPL_TIMING_RPFACT);
#endif
        /*
 * End of HPLAI_pafact
 */
    }

#ifdef __cplusplus
}
#endif
