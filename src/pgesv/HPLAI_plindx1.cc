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
void HPLAI_plindx1(
    HPLAI_T_panel *PANEL,
    const int K,
    const int *IPID,
    int *IPA,
    int *LINDXA,
    int *LINDXAU,
    int *IPLEN,
    int *IPMAP,
    int *IPMAPM1,
    int *PERMU,
    int *IWORK)
#else
void HPLAI_plindx1(PANEL, K, IPID, IPA, LINDXA, LINDXAU, IPLEN, IPMAP, IPMAPM1, PERMU, IWORK)
    HPLAI_T_panel *PANEL;
const int K;
const int *IPID;
int *IPA;
int *LINDXA;
int *LINDXAU;
int *IPLEN;
int *IPMAP;
int *IPMAPM1;
int *PERMU;
int *IWORK;
#endif
{
    /* 
 * Purpose
 * =======
 *
 * HPLAI_plindx1 computes two local arrays  LINDXA and  LINDXAU  containing
 * the  local  source and final destination position  resulting from the
 * application of row interchanges.  In addition, this function computes
 * three arrays IPLEN, IPMAP and IPMAPM1  that contain  the  logarithmic
 * mapping information for the spreading phase.
 *
 * Arguments
 * =========
 *
 * PANEL   (local input/output)          HPLAI_T_panel *
 *         On entry,  PANEL  points to the data structure containing the
 *         panel information.
 *
 * K       (global input)                const int
 *         On entry, K specifies the number of entries in IPID.  K is at
 *         least 2*N, and at most 4*N.
 *
 * IPID    (global input)                const int *
 *         On entry,  IPID  is an array of length K. The first K entries
 *         of that array contain the src and final destination resulting
 *         from the application of the interchanges.
 *
 * IPA     (global output)               int *
 *         On exit,  IPA  specifies  the number of rows that the current
 *         process row has that either belong to U  or should be swapped
 *         with remote rows of A.
 *
 * LINDXA  (global output)               int *
 *         On entry, LINDXA  is an array of dimension 2*N. On exit, this
 *         array contains the local indexes of the rows of A I have that
 *         should be copied into U.
 *
 * LINDXAU (global output)               int *
 *         On exit, LINDXAU  is an array of dimension 2*N. On exit, this
 *         array contains  the local destination  information encoded as
 *         follows.  If LINDXAU(k) >= 0, row  LINDXA(k)  of A  is  to be
 *         copied in U at position LINDXAU(k).  Otherwise, row LINDXA(k)
 *         of A should be locally copied into A(-LINDXAU(k),:).
 *
 * IPLEN   (global output)               int *
 *         On entry, IPLEN is an array of dimension NPROW + 1. On  exit,
 *         this array is such that  IPLEN[i]  is the number of rows of A
 *         in  the  processes  before  process  IPMAP[i]  after the sort
 *         with the convention that IPLEN[nprow]  is the total number of
 *         rows of the panel.  In other words IPLEN[i+1]-IPLEN[i] is the
 *         local number of rows of A that should be moved to the process
 *         IPMAP[i]. IPLEN is such that the number of rows of the source
 *         process  row can be computed as  IPLEN[1] - IPLEN[0], and the
 *         remaining  entries  of  this  array  are  sorted  so that the
 *         quantities IPLEN[i+1] - IPLEN[i] are logarithmically sorted.
 *
 * IPMAP   (global output)               int *
 *         On entry, IPMAP is an array of dimension NPROW. On exit, this
 *         array contains  the logarithmic mapping of the processes.  In
 *         other words, IPMAP[myrow] is the corresponding sorted process
 *         coordinate.
 *
 * IPMAPM1 (global output)               int *
 *         On entry, IPMAPM1  is an array of dimension NPROW.  On  exit,
 *         this  array  contains  the inverse of the logarithmic mapping
 *         contained  in  IPMAP:  IPMAPM1[ IPMAP[i] ] = i,  for all i in
 *         [0.. NPROCS)
 *
 * PERMU   (global output)               int *
 *         On entry,  PERMU  is an array of dimension JB. On exit, PERMU
 *         contains  a sequence of permutations,  that should be applied
 *         in increasing order to permute in place the row panel U.
 *
 * IWORK   (workspace)                   int *
 *         On entry, IWORK is a workarray of dimension 2*JB.
 *
 * ---------------------------------------------------------------------
 */
    /*
 * .. Local Variables ..
 */
    int *iwork;
    int dst, dstrow, fndd, i, ia, icurrow, il,
        ip, ipU, iroff, j, jb, myrow, nb, nprow,
        src, srcrow;
    /* ..
 * .. Executable Statements ..
 */
    /*
 * Logarithmic sort of the processes - compute IPMAP, IPLEN and IPMAPM1
 */
    HPLAI_plindx10(PANEL, K, IPID, IPLEN, IPMAP, IPMAPM1);
    /*
 * Compute the local arrays  LINDXA  and  LINDXAU  containing  the local
 * source and final destination position resulting from  the application
 * of N interchanges. Compute LINDXA and LINDXAU in icurrow,  and LINDXA
 * elsewhere and PERMU in every process.
 */
    myrow = PANEL->grid->myrow;
    nprow = PANEL->grid->nprow;
    jb = PANEL->jb;
    nb = PANEL->nb;
    ia = PANEL->ia;
    iroff = PANEL->ii;
    icurrow = PANEL->prow;

    iwork = IWORK + jb;

    if (myrow == icurrow)
    {
        for (i = 0, ip = 0, ipU = 0; i < K; i += 2)
        {
            src = IPID[i];
            Mindxg2p(src, nb, nb, srcrow, 0, nprow);

            if (srcrow == icurrow)
            {
                dst = IPID[i + 1];
                Mindxg2p(dst, nb, nb, dstrow, 0, nprow);

                Mindxg2l(il, src, nb, nb, myrow, 0, nprow);
                LINDXA[ip] = il - iroff;

                if ((dstrow == icurrow) && (dst - ia < jb))
                {
                    PERMU[ipU] = dst - ia;
                    il = IPMAPM1[dstrow];
                    j = IPLEN[il];
                    iwork[ipU] = LINDXAU[ip] = j;
                    IPLEN[il]++;
                    ipU++;
                }
                else if (dstrow != icurrow)
                {
                    j = 0;
                    do
                    {
                        fndd = (dst == IPID[j]);
                        j += 2;
                    } while (!fndd && (j < K));

                    PERMU[ipU] = IPID[j - 1] - ia;
                    il = IPMAPM1[dstrow];
                    j = IPLEN[il];
                    iwork[ipU] = LINDXAU[ip] = j;
                    IPLEN[il]++;
                    ipU++;
                }
                else if ((dstrow == icurrow) && (dst - ia >= jb))
                {
                    Mindxg2l(il, dst, nb, nb, myrow, 0, nprow);
                    LINDXAU[ip] = iroff - il;
                }
                ip++;
            }
        }
        *IPA = ip;
    }
    else
    {
        for (i = 0, ip = 0, ipU = 0; i < K; i += 2)
        {
            src = IPID[i];
            Mindxg2p(src, nb, nb, srcrow, 0, nprow);
            dst = IPID[i + 1];
            Mindxg2p(dst, nb, nb, dstrow, 0, nprow);
            /*
 * LINDXA[i] is the local index of the row of A that belongs into U
 */
            if (myrow == dstrow)
            {
                Mindxg2l(il, dst, nb, nb, myrow, 0, nprow);
                LINDXA[ip] = il - iroff;
                ip++;
            }
            /*
 * iwork[i] is the local (current) position  index in U
 * PERMU[i] is the local (final) destination index in U
 */
            if (srcrow == icurrow)
            {
                if ((dstrow == icurrow) && (dst - ia < jb))
                {
                    PERMU[ipU] = dst - ia;
                    il = IPMAPM1[dstrow];
                    iwork[ipU] = IPLEN[il];
                    IPLEN[il]++;
                    ipU++;
                }
                else if (dstrow != icurrow)
                {
                    j = 0;
                    do
                    {
                        fndd = (dst == IPID[j]);
                        j += 2;
                    } while (!fndd && (j < K));
                    PERMU[ipU] = IPID[j - 1] - ia;
                    il = IPMAPM1[dstrow];
                    iwork[ipU] = IPLEN[il];
                    IPLEN[il]++;
                    ipU++;
                }
            }
        }
        *IPA = 0;
    }
    /*
 * Simplify iwork and PERMU, return in PERMU the sequence of permutation
 * that need to be apply to U after it has been broadcast.
 */
    HPLAI_perm(jb, iwork, PERMU, IWORK);
    /*
 * Reset IPLEN to its correct value
 */
    for (i = nprow; i > 0; i--)
        IPLEN[i] = IPLEN[i - 1];
    IPLEN[0] = 0;
    /*
 * End of HPLAI_plindx1
 */
}
