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

#ifdef STDC_HEADERS
void HPLAI_palaswp01T(
    HPLAI_T_panel *PBCST,
    int *IFLAG,
    HPLAI_T_panel *PANEL,
    const int NN)
#else
void HPLAI_palaswp01T(PBCST, IFLAG, PANEL, NN)
    HPLAI_T_panel *PBCST;
int *IFLAG;
HPLAI_T_panel *PANEL;
const int NN;
#endif
{
    /* 
 * Purpose
 * =======
 *
 * HPLAI_palaswp01T applies the  NB  row interchanges to  NN columns of the
 * trailing submatrix and broadcast a column panel.
 *  
 * A "Spread then roll" algorithm performs  the swap :: broadcast  of the
 * row panel U at once,  resulting in a minimal communication volume  and
 * a "very good"  use of the connectivity if available.  With  P  process
 * rows  and  assuming  bi-directional links,  the  running time  of this
 * function can be approximated by:
 *  
 *    (log_2(P)+(P-1)) * lat +   K * NB * LocQ(N) / bdwth
 *  
 * where  NB  is the number of rows of the row panel U,  N is the global
 * number of columns being updated,  lat and bdwth  are the latency  and
 * bandwidth  of  the  network  for  HPLAI_T_AFLOAT  precision real words.  K is
 * a constant in (2,3] that depends on the achieved bandwidth  during  a
 * simultaneous  message exchange  between two processes.  An  empirical
 * optimistic value of K is typically 2.4.
 *
 * Arguments
 * =========
 *
 * PBCST   (local input/output)          HPLAI_T_panel *
 *         On entry,  PBCST  points to the data structure containing the
 *         panel (to be broadcast) information.
 *
 * IFLAG   (local input/output)          int *
 *         On entry, IFLAG  indicates  whether or not  the broadcast has
 *         already been completed.  If not,  probing will occur, and the
 *         outcome will be contained in IFLAG on exit.
 *
 * PANEL   (local input/output)          HPLAI_T_panel *
 *         On entry,  PANEL  points to the data structure containing the
 *         panel information.
 *
 * NN      (local input)                 const int
 *         On entry, NN specifies  the  local  number  of columns of the
 *         trailing  submatrix  to  be swapped and broadcast starting at
 *         the current position. NN must be at least zero.
 *
 * ---------------------------------------------------------------------
 */
    /*
 * .. Local Variables ..
 */
    HPLAI_T_AFLOAT *A, *U;
    int *ipID, *iplen, *ipmap, *ipmapm1,
        *iwork, *lindxA = NULL, *lindxAU,
                *permU;
    static int equil = -1;
    int icurrow, *iflag, *ipA, *ipl, jb, k,
        lda, myrow, n, nprow;
#define LDU n
    /* ..
 * .. Executable Statements ..
 */
    n = PANEL->n;
    n = Mmin(NN, n);
    jb = PANEL->jb;
    /*
 * Quick return if there is nothing to do
 */
    if ((n <= 0) || (jb <= 0))
        return;
#ifdef HPL_DETAILED_TIMING
    HPL_ptimer(HPL_TIMING_LASWP);
#endif
    /*
 * Decide whether equilibration should be performed or not
 */
    if (equil == -1)
        equil = PANEL->algo->equil;
    /*
 * Retrieve parameters from the PANEL data structure
 */
    nprow = PANEL->grid->nprow;
    myrow = PANEL->grid->myrow;
    A = PANEL->A;
    U = PANEL->U;
    iflag = PANEL->IWORK;
    lda = PANEL->lda;
    icurrow = PANEL->prow;
    /*
 * Compute ipID (if not already done for this panel). lindxA and lindxAU
 * are of length at most 2*jb - iplen is of size nprow+1, ipmap, ipmapm1
 * are of size nprow,  permU is of length jb, and  this function needs a 
 * workspace of size max( 2 * jb (plindx1), nprow+1(equil)): 
 * 1(iflag) + 1(ipl) + 1(ipA) + 9*jb + 3*nprow + 1 + MAX(2*jb,nprow+1)
 * i.e. 4 + 9*jb + 3*nprow + max(2*jb, nprow+1);
 */
    k = (int)((unsigned int)(jb) << 1);
    ipl = iflag + 1;
    ipID = ipl + 1;
    ipA = ipID + ((unsigned int)(k) << 1);
    lindxA = ipA + 1;
    lindxAU = lindxA + k;
    iplen = lindxAU + k;
    ipmap = iplen + nprow + 1;
    ipmapm1 = ipmap + nprow;
    permU = ipmapm1 + nprow;
    iwork = permU + jb;

    if (*iflag == -1) /* no index arrays have been computed so far */
    {
        HPLAI_pipid(PANEL, ipl, ipID);
        HPLAI_plindx1(PANEL, *ipl, ipID, ipA, lindxA, lindxAU, iplen,
                      ipmap, ipmapm1, permU, iwork);
        *iflag = 1;
    }
    else if (*iflag == 0) /* HPLAI_palaswp00T called before: reuse ipID */
    {
        HPLAI_plindx1(PANEL, *ipl, ipID, ipA, lindxA, lindxAU, iplen,
                      ipmap, ipmapm1, permU, iwork);
        *iflag = 1;
    }
    else if ((*iflag == 1) && (equil != 0))
    { /* HPLAI_palaswp01T was call before only re-compute IPLEN, IPMAP */
        HPLAI_plindx10(PANEL, *ipl, ipID, iplen, ipmap, ipmapm1);
        *iflag = 1;
    }
    /*
 * Copy into U the rows to be spread (local to icurrow)
 */
    if (myrow == icurrow)
    {
        HPLAI_alaswp01T(*ipA, n, A, lda, U, LDU, lindxA, lindxAU);
    }
    /*
 * Spread U - optionally probe for column panel
 */
    HPLAI_spreadT(PBCST, IFLAG, PANEL, blas::Side::Right, n, U, LDU, 0, iplen,
                  ipmap, ipmapm1);
    /*
 * Local exchange (everywhere but in process row icurrow)
 */
    if (myrow != icurrow)
    {
        k = ipmapm1[myrow];
        HPLAI_alaswp06T(iplen[k + 1] - iplen[k], n, A, lda, Mptr(U, 0, iplen[k], LDU), LDU, lindxA);
    }
    /*
 * Equilibration
 */
    if (equil != 0)
        HPLAI_equil(PBCST, IFLAG, PANEL, blas::Op::Trans, n, U, LDU, iplen, ipmap,
                    ipmapm1, iwork);
    /*
 * Rolling phase
 */
    HPLAI_rollT(PBCST, IFLAG, PANEL, n, U, LDU, iplen, ipmap, ipmapm1);
    /*
 * Permute U in every process row
 */
    HPLAI_alaswp10N(n, jb, U, LDU, permU);

#ifdef HPL_DETAILED_TIMING
    HPL_ptimer(HPL_TIMING_LASWP);
#endif
    /*
 * End of HPLAI_palaswp01T
 */
}
