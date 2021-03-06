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
void HPLAI_equil(
    HPLAI_T_panel *PBCST,
    int *IFLAG,
    HPLAI_T_panel *PANEL,
    const blas::Op TRANS,
    const int N,
    HPLAI_T_AFLOAT *U,
    const int LDU,
    int *IPLEN,
    const int *IPMAP,
    const int *IPMAPM1,
    int *IWORK)
#else
void HPLAI_equil(PBCST, IFLAG, PANEL, TRANS, N, U, LDU, IPLEN, IPMAP, IPMAPM1, IWORK)
    HPLAI_T_panel *PBCST;
int *IFLAG;
HPLAI_T_panel *PANEL;
const blas::Op TRANS;
const int N;
HPLAI_T_AFLOAT *U;
const int LDU;
int *IPLEN;
const int *IPMAP;
const int *IPMAPM1;
int *IWORK;
#endif
{
    /* 
 * Purpose
 * =======
 *
 * HPLAI_equil equilibrates  the  local  pieces  of U, so that on exit to
 * this function, pieces of U contained in every process row are of the
 * same size. This phase makes the rolling phase optimal.  In addition,
 * this  function probes  for  the  column panel L and forwards it when
 * possible.
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
 *         panel (to be equilibrated) information.
 *
 * TRANS   (global input)                const blas::Op
 *         On entry, TRANS specifies whether  U  is stored in transposed
 *         or non-transposed form.
 *
 * N       (local input)                 const int
 *         On entry, N  specifies the number of rows or columns of  U. N
 *         must be at least 0.
 *
 * U       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry,  U  is an array of dimension (LDU,*) containing the
 *         local pieces of U in each process row.
 *
 * LDU     (local input)                 const int
 *         On entry, LDU specifies the local leading dimension of U. LDU
 *         should be at least MAX(1,IPLEN[nprow]) when  U  is stored  in
 *         non-transposed form, and MAX(1,N) otherwise.
 *
 * IPLEN   (global input)                int *
 *         On entry, IPLEN is an array of dimension NPROW+1.  This array
 *         is such that IPLEN[i+1] - IPLEN[i] is the number of rows of U
 *         in process IPMAP[i].
 *
 * IPMAP   (global input)                const int *
 *         On entry, IPMAP is an array of dimension  NPROW.  This  array
 *         contains  the  logarithmic mapping of the processes. In other
 *         words, IPMAP[myrow]  is the absolute coordinate of the sorted
 *         process.
 *
 * IPMAPM1 (global input)                const int *
 *         On entry, IPMAPM1  is an array of dimension NPROW. This array
 *         contains  the inverse of the logarithmic mapping contained in
 *         IPMAP: For i in [0.. NPROCS) IPMAPM1[IPMAP[i]] = i.
 *
 * IWORK   (workspace)                   int *
 *         On entry, IWORK is a workarray of dimension NPROW+1.
 *
 * ---------------------------------------------------------------------
 */
    /*
 * .. Local Variables ..
 */
    int i, ip, ipU, ipcur, iprow, iptgt, lastrow,
        left, npm1, nprow, ll, llU, llcur, lltgt,
        right, slen, smax, smin;
    /* ..
 * .. Executable Statements ..
 */
    if ((npm1 = (nprow = PANEL->grid->nprow) - 1) <= 1)
        return;
    /*
 * If the current distribution of the pieces of U is already optimal for
 * the rolling phase, then return imediately.  The  optimal distribution
 * is such that ip processes have smax items and the remaining processes
 * only have smin items. Another way to check this is to verify that all
 * differences IPLEN[i+1] - IPLEN[i] are either smin or smax.
 */
    smax = ((slen = IPLEN[nprow]) + npm1) / nprow;
    ip = slen - nprow * (smin = slen / nprow);

    iprow = 0;
    do
    {
        ll = IPLEN[iprow + 1] - IPLEN[iprow];
        iprow++;
    } while ((iprow < nprow) && ((ll == smin) || (ll == smax)));

    if (iprow == nprow)
        return;
    /*
 * Now,  we are sure  the distribution of the pieces of U is not optimal
 * with respect to the rolling phase,  thus  perform  equilibration.  Go
 * through the list of processes:  Processes  that have rows that do not
 * belong to them  with respect to the optimal mapping spread them  in a
 * logarithmic fashion. To simplify a little bit the implementation, and
 * mainly the packing, a source process row spreads its data to its left
 * first, and then to its right.
 */
    IWORK[nprow] = slen;

    for (iprow = 0; iprow < nprow; iprow++)
    {
        llU = IPLEN[iprow + 1] - (ipU = IPLEN[iprow]);
        if (iprow < ip)
        {
            lltgt = smax;
            iptgt = iprow * smax;
        }
        else
        {
            lltgt = smin;
            iptgt = iprow * smin + ip;
        }

        left = (ipU < iptgt);
        right = (iptgt + lltgt < ipU + llU);
        /*
 * If I have something to spread to either the left or the right
 */
        if ((llU > 0) && (left || right))
        { /* Figure out how much every other process should have */

            ipcur = ipU;
            llcur = llU;

            for (i = 0; i < nprow; i++)
            {
                if (i < ip)
                {
                    lltgt = smax;
                    iptgt = i * smax;
                }
                else
                {
                    lltgt = smin;
                    iptgt = i * smin + ip;
                }
                lastrow = iptgt + lltgt - 1;

                if ((lastrow >= ipcur) && (llcur > 0))
                {
                    ll = lastrow - ipcur + 1;
                    ll = Mmin(ll, llcur);
                    llcur -= ll;
                }
                else
                {
                    ll = 0;
                }

                IWORK[i] = ipcur;
                ipcur += ll;
                IWORK[i + 1] = ipcur;
            }
            /*
 * Equilibration phase
 */
            if (TRANS == blas::Op::NoTrans)
            {
                if (left)
                {
                    HPLAI_spreadN(PBCST, IFLAG, PANEL, blas::Side::Left, N, U, LDU,
                                  iprow, IWORK, IPMAP, IPMAPM1);
                }

                if (right)
                {
                    HPLAI_spreadN(PBCST, IFLAG, PANEL, blas::Side::Right, N, U, LDU,
                                  iprow, IWORK, IPMAP, IPMAPM1);
                }
            }
            else
            {
                if (left)
                {
                    HPLAI_spreadT(PBCST, IFLAG, PANEL, blas::Side::Left, N, U, LDU,
                                  iprow, IWORK, IPMAP, IPMAPM1);
                }

                if (right)
                {
                    HPLAI_spreadT(PBCST, IFLAG, PANEL, blas::Side::Right, N, U, LDU,
                                  iprow, IWORK, IPMAP, IPMAPM1);
                }
            }
        }
    }
    /*
 * Finally update  IPLEN  with the indexes corresponding to the new dis-
 * tribution of U - IPLEN[nprow] remained unchanged.
 */
    for (i = 0; i < nprow; i++)
        IPLEN[i] = (i < ip ? i * smax : i * smin + ip);
    /*
 * End of HPLAI_equil
 */
}
