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
    HPLAI_T_AFLOAT HPLAI_palange(
        const HPLAI_T_grid *GRID,
        const HPLAI_T_NORM NORM,
        const int M,
        const int N,
        const int NB,
        const HPLAI_T_AFLOAT *A,
        const int LDA)
#else
HPLAI_T_AFLOAT HPLAI_palange(GRID, NORM, M, N, NB, A, LDA)
    const HPLAI_T_grid *GRID;
const HPLAI_T_NORM NORM;
const int M;
const int N;
const int NB;
const HPLAI_T_AFLOAT *A;
const int LDA;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_palange returns  the value of the one norm,  or the infinity norm,
 * or the element of largest absolute value of a distributed matrix A:  
 *  
 *  
 *    max(abs(A(i,j))) when NORM = HPLAI_NORM_A,                          
 *    norm1(A),        when NORM = HPLAI_NORM_1,                          
 *    normI(A),        when NORM = HPLAI_NORM_I,                          
 *  
 * where norm1 denotes the one norm of a matrix (maximum column sum) and
 * normI denotes  the infinity norm of a matrix (maximum row sum).  Note
 * that max(abs(A(i,j))) is not a matrix norm.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 const HPLAI_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * NORM    (global input)                const HPLAI_T_NORM
 *         On entry,  NORM  specifies  the  value to be returned by this
 *         function as described above.
 *
 * M       (global input)                const int
 *         On entry,  M  specifies  the number  of rows of the matrix A.
 *         M must be at least zero.
 *
 * N       (global input)                const int
 *         On entry,  N specifies the number of columns of the matrix A.
 *         N must be at least zero.
 *
 * NB      (global input)                const int
 *         On entry,  NB specifies the blocking factor used to partition
 *         and distribute the matrix. NB must be larger than one.
 *
 * A       (local input)                 const HPLAI_T_AFLOAT *
 *         On entry,  A  points to an array of dimension  (LDA,LocQ(N)),
 *         that contains the local pieces of the distributed matrix A.
 *
 * LDA     (local input)                 const int
 *         On entry, LDA specifies the leading dimension of the array A.
 *         LDA must be at least max(1,LocP(M)).
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        HPLAI_T_AFLOAT s, v0 = HPLAI_rzero, *work = NULL;
        MPI_Comm Acomm, Ccomm, Rcomm;
        int ii, jj, mp, mycol, myrow, npcol, nprow,
            nq;
        /* ..
 * .. Executable Statements ..
 */
        (void)HPLAI_grid_info(GRID, &nprow, &npcol, &myrow, &mycol);
        Rcomm = GRID->row_comm;
        Ccomm = GRID->col_comm;
        Acomm = GRID->all_comm;

        Mnumroc(mp, M, NB, NB, myrow, 0, nprow);
        Mnumroc(nq, N, NB, NB, mycol, 0, npcol);

        if (Mmin(M, N) == 0)
        {
            return (v0);
        }
        else if (NORM == HPLAI_NORM_A)
        {
            /*
 * max( abs( A ) )
 */
            if ((nq > 0) && (mp > 0))
            {
                for (jj = 0; jj < nq; jj++)
                {
                    for (ii = 0; ii < mp; ii++)
                    {
                        v0 = Mmax(v0, Mabs(*A));
                        A++;
                    }
                    A += LDA - mp;
                }
            }
            (void)HPLAI_reduce_AFLOAT((void *)(&v0), 1, HPLAI_max_AFLOAT, 0,
                                      Acomm);
        }
        else if (NORM == HPLAI_NORM_1)
        {
            /*
 * Find norm_1( A ).
 */
            if (nq > 0)
            {
                work = (HPLAI_T_AFLOAT *)malloc((size_t)(nq) * sizeof(HPLAI_T_AFLOAT));
                if (work == NULL)
                {
                    HPLAI_pabort(__LINE__, "HPLAI_palange", "Memory allocation failed");
                }

                for (jj = 0; jj < nq; jj++)
                {
                    s = HPLAI_rzero;
                    for (ii = 0; ii < mp; ii++)
                    {
                        s += Mabs(*A);
                        A++;
                    }
                    work[jj] = s;
                    A += LDA - mp;
                }
                /*
 * Find sum of global matrix columns, store on row 0 of process grid
 */
                (void)HPLAI_reduce_AFLOAT((void *)(work), nq, HPLAI_sum_AFLOAT,
                                          0, Ccomm);
                /*
 * Find maximum sum of columns for 1-norm
 */
                if (myrow == 0)
                {
                    v0 = work[blas::iamax<HPLAI_T_AFLOAT>(nq, work, 1)];
                    v0 = Mabs(v0);
                }
                if (work)
                    free(work);
            }
            /*
 * Find max in row 0, store result in process (0,0)
 */
            if (myrow == 0)
                (void)HPLAI_reduce_AFLOAT((void *)(&v0), 1, HPLAI_max_AFLOAT, 0,
                                          Rcomm);
        }
        else if (NORM == HPLAI_NORM_I)
        {
            /*
 * Find norm_inf( A )
 */
            if (mp > 0)
            {
                work = (HPLAI_T_AFLOAT *)malloc((size_t)(mp) * sizeof(HPLAI_T_AFLOAT));
                if (work == NULL)
                {
                    HPLAI_pabort(__LINE__, "HPLAI_palange", "Memory allocation failed");
                }

                for (ii = 0; ii < mp; ii++)
                {
                    work[ii] = HPLAI_rzero;
                }

                for (jj = 0; jj < nq; jj++)
                {
                    for (ii = 0; ii < mp; ii++)
                    {
                        work[ii] += Mabs(*A);
                        A++;
                    }
                    A += LDA - mp;
                }
                /*       
 * Find sum of global matrix rows, store on column 0 of process grid
 */
                (void)HPLAI_reduce_AFLOAT((void *)(work), mp, HPLAI_sum_AFLOAT,
                                          0, Rcomm);
                /*       
 * Find maximum sum of rows for inf-norm
 */
                if (mycol == 0)
                {
                    v0 = work[blas::iamax<HPLAI_T_AFLOAT>(mp, work, 1)];
                    v0 = Mabs(v0);
                }
                if (work)
                    free(work);
            }
            /*
 * Find max in column 0, store result in process (0,0)
 */
            if (mycol == 0)
                (void)HPLAI_reduce_AFLOAT((void *)(&v0), 1, HPLAI_max_AFLOAT,
                                          0, Ccomm);
        }
        /*
 * Broadcast answer to every process in the grid
 */
        (void)HPLAI_broadcast_AFLOAT((void *)(&v0), 1, 0, Acomm);

        return (v0);
        /*
 * End of HPLAI_palange
 */
    }

#ifdef __cplusplus
}
#endif
