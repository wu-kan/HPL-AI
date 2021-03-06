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
void HPLAI_paupdateTT(
    HPLAI_T_panel *PBCST,
    int *IFLAG,
    HPLAI_T_panel *PANEL,
    const int NN)
#else
void HPLAI_paupdateTT(PBCST, IFLAG, PANEL, NN)
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
 * HPLAI_paupdateTT broadcast - forward the panel PBCST and simultaneously
 * applies the row interchanges and updates part of the trailing  (using
 * the panel PANEL) submatrix.
 *
 * Arguments
 * =========
 *
 * PBCST   (local input/output)          HPLAI_T_panel *
 *         On entry,  PBCST  points to the data structure containing the
 *         panel (to be broadcast) information.
 *
 * IFLAG   (local output)                int *
 *         On exit,  IFLAG  indicates  whether or not  the broadcast has
 *         been completed when PBCST is not NULL on entry. In that case,
 *         IFLAG is left unchanged.
 *
 * PANEL   (local input/output)          HPLAI_T_panel *
 *         On entry,  PANEL  points to the data structure containing the
 *         panel (to be updated) information.
 *
 * NN      (local input)                 const int
 *         On entry, NN specifies  the  local  number  of columns of the
 *         trailing  submatrix  to be updated  starting  at the  current
 *         position. NN must be at least zero.
 *
 * ---------------------------------------------------------------------
 */
    /*
 * .. Local Variables ..
 */
    HPLAI_T_AFLOAT *Aptr, *L1ptr, *L2ptr, *Uptr, *dpiv;
    int *ipiv;
#ifdef HPL_CALL_VSIPL
    vsip_mview_d *Av0, *Av1, *Lv0, *Lv1, *Uv0, *Uv1;
#endif
    int curr, i, iroff, jb, lda, ldl2, mp, n, nb,
        nq0, nn, test;
    static int tswap = 0;
    static HPLAI_T_SWAP fswap = HPLAI_NO_SWP;
#define LDU n
/* ..
 * .. Executable Statements ..
 */
#ifdef HPL_DETAILED_TIMING
    HPL_ptimer(HPL_TIMING_UPDATE);
#endif
    nb = PANEL->nb;
    jb = PANEL->jb;
    n = PANEL->nq;
    lda = PANEL->lda;
    if (NN >= 0)
        n = Mmin(NN, n);
    /*
 * There is nothing to update, enforce the panel broadcast.
 */
    if ((n <= 0) || (jb <= 0))
    {
        if (PBCST != NULL)
        {
            do
            {
                (void)HPLAI_bcast(PBCST, IFLAG);
            } while (*IFLAG != HPLAI_SUCCESS);
        }
#ifdef HPL_DETAILED_TIMING
        HPL_ptimer(HPL_TIMING_UPDATE);
#endif
        return;
    }
    /*
 * Enable/disable the column panel probing mechanism
 */
    (void)HPLAI_bcast(PBCST, &test);
    /*
 * 1 x Q case
 */
    if (PANEL->grid->nprow == 1)
    {
        Aptr = PANEL->A;
        L2ptr = PANEL->L2;
        L1ptr = PANEL->L1;
        ldl2 = PANEL->ldl2;
        dpiv = PANEL->DPIV;
        ipiv = PANEL->IWORK;
        mp = PANEL->mp - jb;
        iroff = PANEL->ii;
        nq0 = 0;
#ifdef HPL_CALL_VSIPL
        /*
 * Admit the blocks
 */
        (void)vsip_blockadmit_d(PANEL->Ablock, VSIP_TRUE);
        (void)vsip_blockadmit_d(PANEL->L2block, VSIP_TRUE);
        /*
 * Create the matrix views
 */
        Av0 = vsip_mbind_d(PANEL->Ablock, 0, 1, lda, lda, PANEL->pmat->nq);
        Lv0 = vsip_mbind_d(PANEL->L2block, 0, 1, ldl2, ldl2, jb);
        /*
 * Create the matrix subviews
 */
        Lv1 = vsip_msubview_d(Lv0, 0, 0, mp, jb);
#endif
        for (i = 0; i < jb; i++)
        {
            ipiv[i] = (int)(dpiv[i]) - iroff;
        }
        /*
 * So far we have not updated anything -  test availability of the panel
 * to be forwarded - If detected forward it and finish the update in one
 * step.
 */
        while (test == HPLAI_KEEP_TESTING)
        {
            nn = n - nq0;
            nn = Mmin(nb, nn);
/*
 * Update nb columns at a time
 */
#ifdef HPL_DETAILED_TIMING
            HPL_ptimer(HPL_TIMING_LASWP);
            HPLAI_alaswp00N(jb, nn, Aptr, lda, ipiv);
            HPL_ptimer(HPL_TIMING_LASWP);
#else
            HPLAI_alaswp00N(jb, nn, Aptr, lda, ipiv);
#endif
            blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Side::Left, blas::Uplo::Upper, blas::Op::Trans,
                                                       blas::Diag::Unit, jb, nn, HPLAI_rone, L1ptr, jb, Aptr, lda);
#ifdef HPL_CALL_VSIPL
            /*
 * Create the matrix subviews
 */
            Uv1 = vsip_msubview_d(Av0, PANEL->ii, PANEL->jj + nq0, jb, nn);
            Av1 = vsip_msubview_d(Av0, PANEL->ii + jb, PANEL->jj + nq0, mp, nn);

            vsip_gemp_d(-HPLAI_rone, Lv1, VSIP_MAT_NTRANS, Uv1, VSIP_MAT_NTRANS,
                        HPLAI_rone, Av1);
            /*
 * Destroy the matrix subviews
 */
            (void)vsip_mdestroy_d(Av1);
            (void)vsip_mdestroy_d(Uv1);
#else
            blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, mp, nn,
                                                                       jb, -HPLAI_rone, L2ptr, ldl2, Aptr, lda, HPLAI_rone,
                                                                       Mptr(Aptr, jb, 0, lda), lda);
#endif
            Aptr = Mptr(Aptr, 0, nn, lda);
            nq0 += nn;

            (void)HPLAI_bcast(PBCST, &test);
        }
        /*
 * The panel has been forwarded at that point, finish the update
 */
        if ((nn = n - nq0) > 0)
        {
#ifdef HPL_DETAILED_TIMING
            HPL_ptimer(HPL_TIMING_LASWP);
            HPLAI_alaswp00N(jb, nn, Aptr, lda, ipiv);
            HPL_ptimer(HPL_TIMING_LASWP);
#else
            HPLAI_alaswp00N(jb, nn, Aptr, lda, ipiv);
#endif
            blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Side::Left, blas::Uplo::Upper, blas::Op::Trans,
                                                       blas::Diag::Unit, jb, nn, HPLAI_rone, L1ptr, jb, Aptr, lda);
#ifdef HPL_CALL_VSIPL
            /*
 * Create the matrix subviews
 */
            Uv1 = vsip_msubview_d(Av0, PANEL->ii, PANEL->jj + nq0, jb, nn);
            Av1 = vsip_msubview_d(Av0, PANEL->ii + jb, PANEL->jj + nq0, mp, nn);

            vsip_gemp_d(-HPLAI_rone, Lv1, VSIP_MAT_NTRANS, Uv1, VSIP_MAT_NTRANS,
                        HPLAI_rone, Av1);
            /*
 * Destroy the matrix subviews
 */
            (void)vsip_mdestroy_d(Av1);
            (void)vsip_mdestroy_d(Uv1);
#else
            blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, mp, nn,
                                                                       jb, -HPLAI_rone, L2ptr, ldl2, Aptr, lda, HPLAI_rone,
                                                                       Mptr(Aptr, jb, 0, lda), lda);
#endif
        }
#ifdef HPL_CALL_VSIPL
        /*
 * Destroy the matrix subviews
 */
        (void)vsip_mdestroy_d(Lv1);
        /*
 * Release the blocks
 */
        (void)vsip_blockrelease_d(vsip_mgetblock_d(Lv0), VSIP_TRUE);
        (void)vsip_blockrelease_d(vsip_mgetblock_d(Av0), VSIP_TRUE);
        /*
 * Destroy the matrix views
 */
        (void)vsip_mdestroy_d(Lv0);
        (void)vsip_mdestroy_d(Av0);
#endif
    }
    else /* nprow > 1 ... */
    {
        /*
 * Selection of the swapping algorithm - swap:broadcast U.
 */
        if (fswap == HPLAI_NO_SWP)
        {
            fswap = PANEL->algo->fswap;
            tswap = PANEL->algo->fsthr;
        }

        if ((fswap == HPLAI_SWAP01) ||
            ((fswap == HPLAI_SW_MIX) && (n > tswap)))
        {
            HPLAI_palaswp01T(PBCST, &test, PANEL, n);
        }
        else
        {
            HPLAI_palaswp00T(PBCST, &test, PANEL, n);
        }
        /*
 * Compute redundantly row block of U and update trailing submatrix
 */
        nq0 = 0;
        curr = (PANEL->grid->myrow == PANEL->prow ? 1 : 0);
        Aptr = PANEL->A;
        L2ptr = PANEL->L2;
        L1ptr = PANEL->L1;
        Uptr = PANEL->U;
        ldl2 = PANEL->ldl2;
        mp = PANEL->mp - (curr != 0 ? jb : 0);
#ifdef HPL_CALL_VSIPL
        /*
 * Admit the blocks
 */
        (void)vsip_blockadmit_d(PANEL->Ablock, VSIP_TRUE);
        (void)vsip_blockadmit_d(PANEL->L2block, VSIP_TRUE);
        (void)vsip_blockadmit_d(PANEL->Ublock, VSIP_TRUE);
        /*
 * Create the matrix views
 */
        Av0 = vsip_mbind_d(PANEL->Ablock, 0, 1, lda, lda, PANEL->pmat->nq);
        Lv0 = vsip_mbind_d(PANEL->L2block, 0, 1, ldl2, ldl2, jb);
        Uv0 = vsip_mbind_d(PANEL->Ublock, 0, 1, LDU, LDU, jb);
        /*
 * Create the matrix subviews
 */
        Lv1 = vsip_msubview_d(Lv0, 0, 0, mp, jb);
#endif
        /*
 * Broadcast has not occured yet, spliting the computational part
 */
        while (test == HPLAI_KEEP_TESTING)
        {
            nn = n - nq0;
            nn = Mmin(nb, nn);

            blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Side::Right, blas::Uplo::Upper, blas::Op::NoTrans,
                                                       blas::Diag::Unit, nn, jb, HPLAI_rone, L1ptr, jb, Uptr, LDU);

            if (curr != 0)
            {
#ifdef HPL_CALL_VSIPL
                /*
 * Create the matrix subviews
 */
                Uv1 = vsip_msubview_d(Uv0, nq0, 0, nn, jb);
                Av1 = vsip_msubview_d(Av0, PANEL->ii + jb, PANEL->jj + nq0, mp, nn);

                vsip_gemp_d(-HPLAI_rone, Lv1, VSIP_MAT_NTRANS, Uv1, VSIP_MAT_TRANS,
                            HPLAI_rone, Av1);
                /*
 * Destroy the matrix subviews
 */
                (void)vsip_mdestroy_d(Av1);
                (void)vsip_mdestroy_d(Uv1);
#else
                blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans, mp, nn,
                                                                           jb, -HPLAI_rone, L2ptr, ldl2, Uptr, LDU, HPLAI_rone,
                                                                           Mptr(Aptr, jb, 0, lda), lda);
#endif
                HPLAI_alatcpy(jb, nn, Uptr, LDU, Aptr, lda);
            }
            else
            {
#ifdef HPL_CALL_VSIPL
                /*
 * Create the matrix subviews
 */
                Uv1 = vsip_msubview_d(Uv0, nq0, 0, nn, jb);
                Av1 = vsip_msubview_d(Av0, PANEL->ii, PANEL->jj + nq0, mp, nn);

                vsip_gemp_d(-HPLAI_rone, Lv1, VSIP_MAT_NTRANS, Uv1, VSIP_MAT_TRANS,
                            HPLAI_rone, Av1);
                /*
 * Destroy the matrix subviews
 */
                (void)vsip_mdestroy_d(Av1);
                (void)vsip_mdestroy_d(Uv1);
#else
                blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans, mp, nn,
                                                                           jb, -HPLAI_rone, L2ptr, ldl2, Uptr, LDU, HPLAI_rone,
                                                                           Aptr, lda);
#endif
            }
            Uptr = Mptr(Uptr, nn, 0, LDU);
            Aptr = Mptr(Aptr, 0, nn, lda);
            nq0 += nn;

            (void)HPLAI_bcast(PBCST, &test);
        }
        /*
 * The panel has been forwarded at that point, finish the update
 */
        if ((nn = n - nq0) > 0)
        {
            blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Side::Right, blas::Uplo::Upper, blas::Op::NoTrans,
                                                       blas::Diag::Unit, nn, jb, HPLAI_rone, L1ptr, jb, Uptr, LDU);

            if (curr != 0)
            {
#ifdef HPL_CALL_VSIPL
                /*
 * Create the matrix subviews
 */
                Uv1 = vsip_msubview_d(Uv0, nq0, 0, nn, jb);
                Av1 = vsip_msubview_d(Av0, PANEL->ii + jb, PANEL->jj + nq0, mp, nn);

                vsip_gemp_d(-HPLAI_rone, Lv1, VSIP_MAT_NTRANS, Uv1, VSIP_MAT_TRANS,
                            HPLAI_rone, Av1);
                /*
 * Destroy the matrix subviews
 */
                (void)vsip_mdestroy_d(Av1);
                (void)vsip_mdestroy_d(Uv1);
#else
                blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans, mp, nn,
                                                                           jb, -HPLAI_rone, L2ptr, ldl2, Uptr, LDU, HPLAI_rone,
                                                                           Mptr(Aptr, jb, 0, lda), lda);
#endif
                HPLAI_alatcpy(jb, nn, Uptr, LDU, Aptr, lda);
            }
            else
            {
#ifdef HPL_CALL_VSIPL
                /*
 * Create the matrix subviews
 */
                Uv1 = vsip_msubview_d(Uv0, nq0, 0, nn, jb);
                Av1 = vsip_msubview_d(Av0, PANEL->ii, PANEL->jj + nq0, mp, nn);

                vsip_gemp_d(-HPLAI_rone, Lv1, VSIP_MAT_NTRANS, Uv1, VSIP_MAT_TRANS,
                            HPLAI_rone, Av1);
                /*
 * Destroy the matrix subviews
 */
                (void)vsip_mdestroy_d(Av1);
                (void)vsip_mdestroy_d(Uv1);
#else
                blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans, mp, nn,
                                                                           jb, -HPLAI_rone, L2ptr, ldl2, Uptr, LDU, HPLAI_rone,
                                                                           Aptr, lda);
#endif
            }
        }
#ifdef HPL_CALL_VSIPL
        /*
 * Destroy the matrix subviews
 */
        (void)vsip_mdestroy_d(Lv1);
        /*
 * Release the blocks
 */
        (void)vsip_blockrelease_d(vsip_mgetblock_d(Uv0), VSIP_TRUE);
        (void)vsip_blockrelease_d(vsip_mgetblock_d(Lv0), VSIP_TRUE);
        (void)vsip_blockrelease_d(vsip_mgetblock_d(Av0), VSIP_TRUE);
        /*
 * Destroy the matrix views
 */
        (void)vsip_mdestroy_d(Uv0);
        (void)vsip_mdestroy_d(Lv0);
        (void)vsip_mdestroy_d(Av0);
#endif
    }

    PANEL->A = Mptr(PANEL->A, 0, n, lda);
    PANEL->nq -= n;
    PANEL->jj += n;
    /*
 * return the outcome of the probe  (should always be  HPLAI_SUCCESS,  the
 * panel broadcast is enforced in that routine).
 */
    if (PBCST != NULL)
        *IFLAG = test;
#ifdef HPL_DETAILED_TIMING
    HPL_ptimer(HPL_TIMING_UPDATE);
#endif
    /*
 * End of HPLAI_paupdateTT
 */
}
