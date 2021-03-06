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
    void HPLAI_pagesv0(
        HPLAI_T_grid *GRID,
        HPLAI_T_palg *ALGO,
        HPLAI_T_pmat *A)
#else
void HPLAI_pagesv0(GRID, ALGO, A)
    HPLAI_T_grid *GRID;
HPLAI_T_palg *ALGO;
HPLAI_T_pmat *A;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_pagesv0 factors a N+1-by-N matrix using LU factorization with row
 * partial pivoting.  The main algorithm  is the "right looking" variant
 * without look-ahead. The lower triangular factor is left unpivoted and
 * the pivots are not returned. The right hand side is the N+1 column of
 * the coefficient matrix.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPLAI_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * ALGO    (global input)                HPLAI_T_palg *
 *         On entry,  ALGO  points to  the data structure containing the
 *         algorithmic parameters.
 *
 * A       (local input/output)          HPLAI_T_pmat *
 *         On entry, A points to the data structure containing the local
 *         array information.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        HPLAI_T_panel **panel = NULL;
        HPLAI_T_UPD_FUN HPLAI_paupdate;
        int N, j, jb, n, nb, tag = MSGID_BEGIN_FACT,
                             test = HPL_KEEP_TESTING;
#ifdef HPL_PROGRESS_REPORT
        double start_time, time, gflops;
#endif
        /* ..
 * .. Executable Statements ..
 */
        if ((N = A->n) <= 0)
            return;

#ifdef HPL_PROGRESS_REPORT
        start_time = HPL_timer_walltime();
#endif

        HPLAI_paupdate = ALGO->upfun;
        nb = A->nb;
        /*
 * Allocate a panel list of length 1 - Allocate panel[0] resources
 */
        panel = (HPLAI_T_panel **)malloc(sizeof(HPLAI_T_panel *));
        if (panel == NULL)
        {
            HPL_pabort(__LINE__, "HPLAI_pagesv0", "Memory allocation failed");
        }

        HPLAI_papanel_new(GRID, ALGO, N, N + 1, Mmin(N, nb), A, 0, 0, tag,
                          &panel[0]);
        /*
 * Loop over the columns of A
 */
        for (j = 0; j < N; j += nb)
        {
            n = N - j;
            jb = Mmin(n, nb);
#ifdef HPL_PROGRESS_REPORT
            /* if this is process 0,0 and not the first panel */
            if (GRID->myrow == 0 && GRID->mycol == 0 && j > 0)
            {
                time = HPL_timer_walltime() - start_time;
                gflops = 2.0 * (N * (double)N * N - n * (double)n * n) / 3.0 / (time > 0.0 ? time : 1e-6) / 1e9;
                HPL_fprintf(stdout, "Column=%09d Fraction=%4.1f%% Gflops=%9.3e\n", j, j * 100.0 / N, gflops);
            }
#endif
            /*
 * Release panel resources - re-initialize panel data structure
 */
            (void)HPLAI_papanel_free(panel[0]);
            HPLAI_papanel_init(GRID, ALGO, n, n + 1, jb, A, j, j, tag, panel[0]);
            /*
 * Factor and broadcast current panel - update
 */
            HPLAI_pafact(panel[0]);
            (void)HPLAI_binit(panel[0]);
            do
            {
                (void)HPLAI_bcast(panel[0], &test);
            } while (test != HPL_SUCCESS);
            (void)HPLAI_bwait(panel[0]);
            HPLAI_paupdate(NULL, NULL, panel[0], -1);
            /*
 * Update message id for next factorization
 */
            tag = MNxtMgid(tag, MSGID_BEGIN_FACT, MSGID_END_FACT);
        }
        /*
 * Release panel resources and panel list
 */
        (void)HPLAI_papanel_disp(&panel[0]);

        if (panel)
            free(panel);
        /*
 * End of HPLAI_pagesv0
 */
    }

#ifdef __cplusplus
}
#endif
