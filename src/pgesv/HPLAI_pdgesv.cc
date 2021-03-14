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

// https://github.com/schuangs/hpl-ai-with-IR/blob/master/src/pgesv/HPL_pLdtrsv.c

static void HPL_pLdtrsv(
    HPL_T_grid *GRID,
    HPL_T_pmat *AMAT)
{
    /* 
 * Purpose
 * =======
 *
 * HPL_pdtrsv solves an lower triangular system of linear equations.
 *  
 * The rhs is the last column of the N by N+1 matrix A. The solve starts
 * in the process  column owning the  1th  column of A, so the rhs b may
 * need to be moved one process column to the left at the beginning. The
 * routine therefore needs  a column  vector in every process column but
 * the one owning  b. The result is  replicated in all process rows, and
 * returned in XR, i.e. XR is of size nq = LOCq( N ) in all processes.
 *  
 * The algorithm uses decreasing one-ring broadcast in process rows  and
 * columns  implemented  in terms of  synchronous communication point to
 * point primitives.  The  lookahead of depth 1 is used to minimize  the
 * critical path. This entire operation is essentially ``latency'' bound
 * and an estimate of its running time is given by:
 *  
 *    (move rhs) lat + N / ( P bdwth ) +            
 *    (solve)    ((N / NB)-1) 2 (lat + NB / bdwth) +
 *               gam2 N^2 / ( P Q ),                
 *  
 * where  gam2   is an estimate of the   Level 2 BLAS rate of execution.
 * There are  N / NB  diagonal blocks. One must exchange  2  messages of
 * length NB to compute the next  NB  entries of the vector solution, as
 * well as performing a total of N^2 floating point operations.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPL_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * AMAT    (local input/output)          HPL_T_pmat *
 *         On entry,  AMAT  points  to the data structure containing the
 *         local array information.
 *
 * ---------------------------------------------------------------------
 */
    /*
 * .. Local Variables ..
 */
    MPI_Comm Ccomm, Rcomm;
    double *A = NULL, *Aprev = NULL, *Aptr, *XC = NULL,
           *XR = NULL, *Xd = NULL, *Xdprev = NULL,
           *W = NULL;
    int Alcol, Alrow, Anpprev, Anp, Anq, Bcol,
        Cmsgid, GridIsNotPx1, GridIsNot1xQ, Rmsgid,
        Wfr = 0, colprev, kb, kbprev, lda, mycol,
        myrow, n, n1, n1p, n1pprev = 0, nb, npcol,
        nprow, rowprev, tmp1, tmp2, np, nq, N;
/* ..
 * .. Executable Statements ..
 */
#ifdef HPL_DETAILED_TIMING
    HPL_ptimer(HPL_TIMING_PTRSV);
#endif
    if ((N = n = AMAT->n) <= 0)
        return;
    nb = AMAT->nb;
    lda = AMAT->ld;
    A = AMAT->A;
    XR = AMAT->X;

    (void)HPL_grid_info(GRID, &nprow, &npcol, &myrow, &mycol);
    Rcomm = GRID->row_comm;
    Rmsgid = MSGID_BEGIN_PTRSV;
    Ccomm = GRID->col_comm;
    Cmsgid = MSGID_BEGIN_PTRSV + 1;
    GridIsNot1xQ = (nprow > 1);
    GridIsNotPx1 = (npcol > 1);
    /*
 * Move the rhs in the process column owning the first column of A.
 */
    /* np and nq are the local dimensions of of A */
    Mnumroc(np, n, nb, nb, myrow, 0, nprow);
    Mnumroc(nq, n, nb, nb, mycol, 0, npcol);

    Anp = 0;
    Anq = 0;

    Alrow = 0;
    Alcol = 0;
    kb = Mmin(n, nb);

    Aptr = (double *)(A);
    XC = Mptr(Aptr, 0, nq, lda);
    Mindxg2p(n, nb, nb, Bcol, 0, npcol);

    if ((np > 0) && (Alcol != Bcol))
    {
        if (mycol == Bcol)
        {
            (void)HPL_send(XC, np, Alcol, Rmsgid, Rcomm);
        }
        else if (mycol == Alcol)
        {
            (void)HPL_recv(XC, np, Bcol, Rmsgid, Rcomm);
        }
    }
    Rmsgid = (Rmsgid + 2 >
                      MSGID_END_PTRSV
                  ? MSGID_BEGIN_PTRSV
                  : Rmsgid + 2);

    /* other rhs area set to 0 */
    if (mycol != Alcol)
    {
        for (tmp1 = 0; tmp1 < np; tmp1++)
            XC[tmp1] = HPL_rzero;
    }
    /*
 * Set up lookahead, that is to perform iteration on the first block column
 */

    /* n1 is the global distance between two neighboring blocks of the same
       process in a block row */
    n1 = (npcol - 1) * nb;
    n1 = Mmax(n1, nb);

    /* allocate work space */
    if (np > 0)
    {
        W = (double *)malloc((size_t)(Mmin(n1, np)) * sizeof(double));
        if (W == NULL)
        {
            HPLAI_pabort(__LINE__, "HPL_pLdtrsv", "Memory allocation failed");
        }
        Wfr = 1;
    }

    /* previous version of parameters */
    Anpprev = Anp;
    Xdprev = Xd = XR;
    Aprev = Aptr;
    tmp1 = kb;
    tmp2 = Mmin(n - kb, n1);
    MnumrocI(n1pprev, tmp2, Mmin(N, tmp1), nb, nb, myrow, 0, nprow);

    /* update local variables */
    if (myrow == Alrow)
    {
        Anp += kb;
    }
    if (mycol == Alcol)
    {
        if (myrow == Alrow)
        {
            blas::trsv<double, double>(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::NoTrans, blas::Diag::Unit,
                                       kb, Aptr, lda, XC, 1);
            blas::copy<double, double>(kb, XC, 1, Xd, 1);
        }

        /* update local variables */
        Xd += kb;
        Anq += kb;
        Aptr += lda * kb;
    }

    /* update global variables */
    rowprev = Alrow;
    Alrow = MModAdd1(Alrow, nprow);
    colprev = Alcol;
    Alcol = MModAdd1(Alcol, npcol);
    kbprev = kb;
    n -= kb;
    tmp1 = N - n + (kb = Mmin(n, nb));
    tmp2 = Mmin(n - kb, n1);
    MnumrocI(n1p, tmp2, Mmin(N, tmp1), nb, nb, myrow, 0, nprow);
    /*
 * Start the operations
 */
    while (n > 0)
    {
        /*
 * Broadcast  (decreasing-ring)  of  previous solution block in previous
 * process column,  compute  partial update of current block and send it
 * to current process column.
 */
        if (mycol == colprev)
        {
            /*
 * Send previous solution block in process row below
 */
            if (myrow == rowprev)
            {
                if (GridIsNot1xQ)
                    (void)HPL_send(Xdprev, kbprev, MModAdd1(myrow, nprow),
                                   Cmsgid, Ccomm);
            }
            else
            {
                (void)HPL_recv(Xdprev, kbprev, MModSub1(myrow, nprow),
                               Cmsgid, Ccomm);
            }
            /*
 * Compute partial update of previous solution block and send it to cur-
 * rent column
 */
            if (n1pprev > 0)
            {
                blas::gemv<double, double, double>(blas::Layout::ColMajor, blas::Op::NoTrans, n1pprev, kbprev,
                                                   -HPL_rone, Aprev + Anp, lda, Xdprev, 1, HPL_rone,
                                                   XC + Anp, 1);
                if (GridIsNotPx1)
                    (void)HPL_send(XC + Anp, n1pprev, Alcol, Rmsgid, Rcomm);
            }
            /*
 * Finish  the (decreasing-ring) broadcast of the solution block in pre-
 * vious process column
 */
            if ((myrow != rowprev) &&
                (myrow != MModSub1(rowprev, nprow)))
                (void)HPL_send(Xdprev, kbprev, MModAdd1(myrow, nprow),
                               Cmsgid, Ccomm);
        }
        else if (mycol == Alcol)
        {
            /*
 * Current  column  receives  and accumulates partial update of previous
 * solution block
 */
            if (n1pprev > 0)
            {
                (void)HPL_recv(W, n1pprev, colprev, Rmsgid, Rcomm);
                blas::axpy<double, double>(n1pprev, HPL_rone, W, 1, XC + Anp, 1);
            }
        }
        /*
 * Solve current diagonal block 
 */
        if ((mycol == Alcol) && (myrow == Alrow))
        {
            blas::trsv<double, double>(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::NoTrans, blas::Diag::Unit,
                                       kb, Aptr + Anp, lda, XC + Anp, 1);
            blas::copy<double, double>(kb, XC + Anp, 1, Xd, 1);
        }
        /*
 *  Finish previous update
 */
        if ((mycol == colprev) && ((tmp1 = Anp + n1pprev) < np))
            blas::gemv<double, double, double>(blas::Layout::ColMajor, blas::Op::NoTrans, np - tmp1, kbprev, -HPL_rone,
                                               Aprev + tmp1, lda, Xdprev, 1, HPL_rone, XC + tmp1, 1);

        /*
 *  Save info of current step and update info for the next step
 */
        if (mycol == Alcol)
        {
            Aprev = Aptr;
            Aptr += lda * kb;
            Anq += kb;
            Xdprev = Xd;
            Xd = XR + Anq;
        }
        if (myrow == Alrow)
        {
            Anpprev = Anp;
            Anp += kb;
        }

        rowprev = Alrow;
        colprev = Alcol;
        n1pprev = n1p;
        kbprev = kb;
        n -= kb;
        Alrow = MModAdd1(Alrow, nprow);
        Alcol = MModAdd1(Alcol, npcol);

        tmp1 = N - n + (kb = Mmin(n, nb));
        tmp2 = Mmin(n - kb, n1);
        MnumrocI(n1p, tmp2, Mmin(N, tmp1), nb, nb, myrow, 0, nprow);

        Rmsgid = (Rmsgid + 2 > MSGID_END_PTRSV ? MSGID_BEGIN_PTRSV : Rmsgid + 2);
        Cmsgid = (Cmsgid + 2 > MSGID_END_PTRSV ? MSGID_BEGIN_PTRSV + 1 : Cmsgid + 2);
    }
    /*
 * Replicate last solution block
 */
    if (mycol == colprev)
        (void)HPL_broadcast((void *)Xdprev, kbprev, HPL_DOUBLE, rowprev,
                            Ccomm);

    if (Wfr)
        free(W);
#ifdef HPL_DETAILED_TIMING
    HPL_ptimer(HPL_TIMING_PTRSV);
#endif
    /*
 * End of HPL_pdtrsv
 */
}

// https://github.com/schuangs/hpl-ai-with-IR/blob/master/src/pir/HPL_pgmres.c

/*
 *  by Junkang Huang, Dec. 2020
 * 
 *  based on the implementation of parallel GMRES:
 *      Parallelization Of The GMRES ---by Morgan Görtz, Lund University.
 *  
 *  and the pioneer work implementing Householder Transformations into GMRES:
 *      Implementation Of The GMRES Method Using Householder Transformations Method
 *      ---by Homer F. Walker, 1988
 */

static double sign(double x)
{
    return x < 0 ? -1 : 1;
}

/* 
 * givens_rotation():
 * 
 * 1. perform:
 *   v <- Jk-1Jk-2...J1J0v
 * 
 * 2. solve for Jk:
 *   s.t. Jkv = (v0,v1,...,vk-1,somevalue,0,...,0)
 * 
 * 3. perform:
 *   v <- Jkv
 *   w <- Jkw
 * 
 * 4. append v to R, that is:
 *   R = [R, v]
 */
static void givens_rotations(
    HPL_T_grid *GRID, /* processes grid information */
    HPL_T_pmat *A,    /* local A */
    double *v,        /* kth column of H */
    double *w,        /* rhs */
    double *R,        /* R matrix */
    double *sinus,    /* sin(theta) */
    double *cosus,    /* cos(theta) */
    const int k,      /* offset */
    const int MM      /* restart size */
)
{
    /* local variables */
    int pi, pi1, ii, ii1, i;
    double tmp;

    /* update v */
    for (i = 0; i < k; ++i)
    {
        /* calculate the process row which possess v[i] and v[i+1] */
        HPL_indxg2lp(&ii, &pi, i, A->nb, A->nb, 0, GRID->nprow);
        HPL_indxg2lp(&ii1, &pi1, i + 1, A->nb, A->nb, 0, GRID->nprow);

        if (pi == pi1)
        { /* if two elememts in one process row, just perform local update */
            if (GRID->myrow == pi)
            {
                /* update v */
                tmp = cosus[i] * v[ii] - sinus[i] * v[ii1];
                v[ii1] = sinus[i] * v[ii] + cosus[i] * v[ii1];
                v[ii] = tmp;
            }
        }
        else
        { /* if two elememts in different process rows, some communication required */
            if (GRID->myrow == pi)
            {
                /* update v */
                HPL_send(&v[ii], 1, pi1, 0, GRID->col_comm);
                HPL_recv(&tmp, 1, pi1, 0, GRID->col_comm);
                v[ii] = cosus[i] * v[ii] - sinus[i] * tmp;
            }
            if (GRID->myrow == pi1)
            {
                /* update v */
                HPL_recv(&tmp, 1, pi, 0, GRID->col_comm);
                HPL_send(&v[ii1], 1, pi, 0, GRID->col_comm);
                v[ii1] = sinus[i] * tmp + cosus[i] * v[ii1];
            }
        }
    }

    /* solve for Jk */

    /* calculate the process row which possess v[k] and v[k+1] */
    HPL_indxg2lp(&ii, &pi, k, A->nb, A->nb, 0, GRID->nprow);
    HPL_indxg2lp(&ii1, &pi1, k + 1, A->nb, A->nb, 0, GRID->nprow);

    if (pi == pi1)
    { /* if two elements in the same process row */
        if (GRID->myrow == pi)
        {
            if (v[ii] * v[ii] + v[ii1] * v[ii1] == 0)
            {
                printf("Error: divided by zero in givens_rotations()\n");
                return;
            }
            /* calculate sin and cos for Jk */
            cosus[k] = v[ii] / sqrt(v[ii] * v[ii] + v[ii1] * v[ii1]);
            sinus[k] = -v[ii1] / sqrt(v[ii] * v[ii] + v[ii1] * v[ii1]);

            /* update v */
            v[ii] = cosus[k] * v[ii] - sinus[k] * v[ii1];
            v[ii1] = 0;
        }
    }
    else
    { /* if two elements not in the same process row */

        /* calculate sin and cos for Jk */
        if (GRID->myrow == pi)
        {
            HPL_recv(&tmp, 1, pi1, 1, GRID->col_comm);
            if (v[ii] * v[ii] + tmp * tmp == 0)
            {
                printf("Error: divided by zero in givens_rotations()\n");
                return;
            }
            cosus[k] = v[ii] / sqrt(v[ii] * v[ii] + tmp * tmp);
            sinus[k] = -tmp / sqrt(v[ii] * v[ii] + tmp * tmp);
            v[ii] = cosus[k] * v[ii] - sinus[k] * tmp;
        }
        if (GRID->myrow == pi1)
        {
            HPL_send(&v[ii1], 1, pi, 1, GRID->col_comm);
            v[ii1] = 0;
        }
    }

    /* broadcast sin and cos */
    HPL_broadcast(&cosus[k], 1, HPL_DOUBLE, pi, GRID->col_comm);
    HPL_broadcast(&sinus[k], 1, HPL_DOUBLE, pi, GRID->col_comm);

    /* update w */
    tmp = cosus[k] * w[k] - sinus[k] * w[k + 1];
    w[k + 1] = sinus[k] * w[k] + cosus[k] * w[k + 1];
    w[k] = tmp;

    /* update R */
    for (i = 0; i < k + 1; ++i)
    {
        HPL_indxg2lp(&ii, &pi, i, A->nb, A->nb, 0, GRID->nprow);
        if (pi == 0)
        { /* if v[i] already in process row 0 */
            if (GRID->myrow == 0)
            {
                /* just perform local R update on process row 0 */
                *Mptr(R, i, k, MM) = v[ii];
            }
        }
        else
        { /* if v[i] is in another process row */
            if (GRID->myrow == pi)
            {
                /* send v[i] to process row 0, i+3 is just a tag 
                    in case of message mismatch */
                HPL_send(&v[ii], 1, 0, i + 3, GRID->col_comm);
            }
            if (GRID->myrow == 0)
            {
                /* process row 0 receive v[i], and update local R */
                HPL_recv(Mptr(R, i, k, MM), 1, pi, i + 3, GRID->col_comm);
            }
        }
    }
    /* broadcast R in process row 0 to all */
    HPL_broadcast(Mptr(R, 0, k, MM), k + 1, HPL_DOUBLE, 0, GRID->col_comm);

    /* end of givens_rotations() */
}

/* 
 * generateHouseholder():
 *  
 * solve for Householder vector u:
 *   s.t. Pkx = (I-2uuT)x = [x0,x1,..,xk-1,alpha,0,..0], alpha != 0
 */
static void generateHouseholder(
    HPL_T_grid *GRID, /* processes grid information */
    HPL_T_pmat *A,    /* local A */
    const double *x,  /* local object vector pointer */
    double *u,        /* local result Householder Vector */
    const int k,      /* order of the Householder */
    double *alpha     /* result variable */
)
{
    /* local variables */
    const int myrow = GRID->myrow;
    int i, ig, mp = A->mp, pi;
    double r = 0;

    for (i = 0; i < mp; ++i)
    {
        ig = HPL_indxl2g(i, A->nb, A->nb, GRID->myrow, 0, GRID->nprow);
        if (ig >= k)
        {
            /* load u[k:] with x[k:] for further operation */
            u[i] = x[i];
            /* calculate (xk*xk) + (xk+1*xk+1) + ...*/
            r += x[i] * x[i];
        }
        else
        {
            /* u[:k] should be 0 */
            u[i] = 0;
        }
    }
    HPL_indxg2lp(&i, &pi, k, A->nb, A->nb, 0, GRID->nprow);
    /* Get the total r on process row which possess u[k] */
    HPL_reduce(&r, 1, HPL_DOUBLE, HPL_sum, pi, GRID->col_comm);

    if (myrow == pi)
    {
        /* perform computation on process pi */
        /* calculate alpha and r */
        *alpha = -sign(x[i]) * sqrt(r);
        r = sqrt(0.5 * ((*alpha) * (*alpha) - x[i] * (*alpha)));
        /* compute the first nonzero value of the transformation vector u */
        u[i] = x[i] - *alpha;
    }
    /* send r and alpha to all process */
    HPL_broadcast(&r, 1, HPL_DOUBLE, pi, GRID->col_comm);
    HPL_broadcast(alpha, 1, HPL_DOUBLE, pi, GRID->col_comm);
    /* apply 1/2r on u for all processes */
    for (i = 0; i < mp; ++i)
    {
        u[i] /= 2. * r;
    }

    /* end of generateHouseholder() */
}

/*
 * applyHouseholder()
 * 
 * perform:
 *    y = Pkx = (I-2uuT)x = x - 2u(x, u)
 */
static void applyHouseholder(
    HPL_T_grid *GRID, /* processes grid information */
    HPL_T_pmat *A,    /* local A */
    const double *x,  /* local object vector pointer */
    const double *u,  /* local result Householder Vector */
    const int k,      /* order of the Householder */
    double *y         /* target vector */
)
{
    /* local variables */
    double segsum = 0;
    int i, ig;

    /* calculate (x, u) */
    for (i = 0; i < A->mp; ++i)
    {
        ig = HPL_indxl2g(i, A->nb, A->nb, GRID->myrow, 0, GRID->nprow);
        if (ig >= k)
        {
            segsum += u[i] * x[i];
        }
    }
    /* sum (x, u) */
    HPL_all_reduce(&segsum, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm);

    /* calculate y = x - 2u(x, u) */
    for (i = 0; i < A->mp; ++i)
    {
        ig = HPL_indxl2g(i, A->nb, A->nb, GRID->myrow, 0, GRID->nprow);
        if (ig >= k)
        {
            y[i] = x[i] - 2 * segsum * u[i];
        }
        else
        {
            y[i] = x[i];
        }
    }

    /* end of applyHouseholder() */
}

/*
 * redB2X()
 * 
 * when performing A*v, if v is distributed along different rows like b, some redistributions
 * needed to perform to make v distributed along different columns like x.
 */
static void redB2X(
    HPL_T_grid *GRID,
    HPL_T_pmat *A,   /* local A */
    const double *v, /* the vector to be redistributed, size: mp */
    double *vc       /* the target space, size: nq */
)
{
    int ig, i, j, jp;
    for (i = 0; i < A->nq - 1; ++i)
    {
        /* find the global index of local vc[i] */
        ig = HPL_indxl2g(i, A->nb, A->nb, GRID->mycol, 0, GRID->npcol);

        /* find the process row and local index of the element vc[i] stored in v */
        HPL_indxg2lp(&j, &jp, ig, A->nb, A->nb, 0, GRID->nprow);

        /* there is one and only one process who contains both vc[i] and v[j] in 
            this process column */
        if (GRID->myrow == jp)
        {
            /* perform local replication */
            vc[i] = v[j];
        }
        /* broadcast the correct vc to other processes in the column */
        HPL_broadcast(&vc[i], 1, HPL_DOUBLE, jp, GRID->col_comm);
    }
}

/*
 * redX2B()
 * 
 * v is distributed along different columns like x, some redistributions
 * needed to perform to make v distributed along different columns like b.
 */
static void redX2B(
    HPL_T_grid *GRID,
    HPL_T_pmat *A,   /* local A */
    const double *v, /* the vector to be redistributed, size: nq */
    double *vc       /* the target space, size: mp */
)
{
    int ig, i, j, jp;
    for (i = 0; i < A->mp; ++i)
    {
        /* find the global index of local vc[i] */
        ig = HPL_indxl2g(i, A->nb, A->nb, GRID->myrow, 0, GRID->nprow);

        /* find the process column and local index of the element vc[i] stored in v */
        HPL_indxg2lp(&j, &jp, ig, A->nb, A->nb, 0, GRID->npcol);

        /* there is one and only one process who contains both vc[i] and v[j] in */
        if (GRID->mycol == jp)
        {
            /* perform local replication */
            vc[i] = v[j];
        }
        /* broadcast the correct vc to other processes in the rows */
        HPL_broadcast(&vc[i], 1, HPL_DOUBLE, jp, GRID->row_comm);
    }
}

/*
 *  HPL_pgmres():
 * 
 */
static int HPL_pgmres(
    HPL_T_grid *GRID,
    HPL_T_pmat *A,       /* local A */
    HPL_T_pmat *factors, /* local LU factors */
    const double *b,     /* local rhs */
    double *x,           /* local solution vector */
    double TOL,          /* tolerance of residual */
    const int MM,        /* restart size */
    const int MAXIT      /* maximum # of total iteration */
)
{
    int prec = 1; /* whether or not to precondition, for debugging */
    /* local variables */
    int i, j, k = 0, start, ready = 0, index, pindex, tarcol;
    double norm, currenterror, tmp;
    int mp = A->mp, nq = A->nq - 1;

    /* bptr point to the rhs area of factors */
    double *bptr = Mptr(factors->A, 0, nq, factors->ld);

    /* distributed storages: each process row stores a part of data */
    double *v = (double *)malloc(mp * sizeof(double));
    double *u = (double *)malloc(mp * sizeof(double));
    double *xt = (double *)malloc(nq * sizeof(double));
    double *H = (double *)malloc(mp * (MM + 1) * sizeof(double));
    double *rhs = (double *)malloc(mp * sizeof(double));

    /* replicated storage: all processes store the whole data */
    double *cosus = (double *)malloc((MM + 1) * sizeof(double));
    double *sinus = (double *)malloc((MM + 1) * sizeof(double));
    double *w = (double *)malloc((MM + 1) * sizeof(double));
    double *R = (double *)malloc(MM * (MM + 1) * sizeof(double));

    /* precondition b into rhs, that is: rhs = U-1L-1b */
    tarcol = HPL_indxg2p(factors->n, factors->nb, factors->nb, 0, GRID->npcol);
    if (prec)
    {
        /* solve Lx = b, then x = L-1b, x is returned at factors->X, which is replicated in process rows */
        if (GRID->mycol == tarcol)
        {
            memcpy(bptr, b, mp * sizeof(double));
        }
        HPL_pLdtrsv(GRID, factors);

        /* redistribute x into column-distributing pattern for next pdtrsv */
        redX2B(GRID, factors, factors->X, rhs);

        /* solve Ux = b, then x = U-1b */
        if (GRID->mycol == tarcol)
        {
            memcpy(bptr, rhs, mp * sizeof(double));
        }
        HPL_pdtrsv(GRID, factors);
        redX2B(GRID, factors, factors->X, rhs);
    }
    else
    {
        memcpy(rhs, b, mp * sizeof(double));
    }

    /* no initial guess so the first residual r0 is just b */
    memcpy(v, rhs, mp * sizeof(double));

    norm = 0;
    /* calculate the norm of r0 */
    for (i = 0; i < mp; ++i)
    {
        norm += v[i] * v[i];
    }
    HPL_all_reduce(&norm, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm);
    norm = sqrt(norm);

    /* store the current error */
    currenterror = norm;

    /* creating w */
    memset(w, 0, (MM + 1) * sizeof(double));
    /* generate and apply the first Housholder transformation
        Householder vector stored in u */
    generateHouseholder(GRID, A, v, u, 0, &w[0]);

    /* stop if approximation is good enough, zero solution returned. */
    if (currenterror < TOL)
    {
        ready = 1;
    }

    /* ------------------------------------------------- */
    /*  Restart Iterations                               */
    /* ------------------------------------------------- */
    for (start = 0; (start < MAXIT) && !ready; ++start)
    {
        /* do the same as above to start the method */
        if (start)
        {
            /* there is initial guess stored in x here from last iteration */
            /* calculate v = Ax */
            blas::gemv<double, double, double>(blas::Layout::ColMajor, blas::Op::NoTrans, mp, nq, HPL_rone,
                                               A->A, A->ld, x, 1, 0, v, 1);
            HPL_all_reduce(v, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);

            /* preconditioning A */
            if (prec)
            {
                if (GRID->mycol == tarcol)
                {
                    memcpy(bptr, v, mp * sizeof(double));
                }
                HPL_pLdtrsv(GRID, factors);
                redX2B(GRID, factors, factors->X, v);

                if (GRID->mycol == tarcol)
                {
                    memcpy(bptr, v, mp * sizeof(double));
                }
                HPL_pdtrsv(GRID, factors);
                redX2B(GRID, factors, factors->X, v);
            }

            /* v = rhs - v = rhs - Ax */
            for (i = 0; i < mp; ++i)
            {
                v[i] = rhs[i] - v[i];
            }
            /* generate P0v = [alpha,0,..,0] */
            memset(w, 0, (MM + 1) * sizeof(double));
            generateHouseholder(GRID, A, v, u, 0, &w[0]);
        }
        /* ------------------------------------------------ */
        /* Householder transformations and Givens rotations */
        /* ------------------------------------------------ */
        for (k = 0; k < MM; ++k)
        {
            if (k >= A->n - 1)
            {
                --k;
                break;
            }
            /* store the current trasformation vector u in H */
            memcpy(Mptr(H, 0, k, mp), u, mp * sizeof(double));

            /* load ek into v */
            memset(v, 0, mp * sizeof(double));
            HPL_indxg2lp(&index, &pindex, k, A->nb, A->nb, 0, GRID->nprow);
            if (GRID->myrow == pindex)
            {
                v[index] = 1;
            }
            /* apply the last k + 1 Householder transformations in reverse order:
                that is : v = P0P1..Pkv */
            for (i = k; i >= 0; --i)
            {
                applyHouseholder(GRID, A, v, Mptr(H, 0, i, mp), i, v);
            }

            /* calculate v = AP0P1..Pkv */
            redB2X(GRID, A, v, xt);
            blas::gemv<double, double, double>(blas::Layout::ColMajor, blas::Op::NoTrans, mp, nq, HPL_rone,
                                               A->A, A->ld, xt, 1, 0, v, 1);
            HPL_all_reduce(v, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);

            /* preconditioning A */
            if (prec)
            {
                if (GRID->mycol == tarcol)
                {
                    memcpy(bptr, v, mp * sizeof(double));
                }
                HPL_pLdtrsv(GRID, factors);
                redX2B(GRID, factors, factors->X, v);

                if (GRID->mycol == tarcol)
                {
                    memcpy(bptr, v, mp * sizeof(double));
                }
                HPL_pdtrsv(GRID, factors);
                redX2B(GRID, factors, factors->X, v);
            }

            /* apply last k + 1 Householder transformations: 
                that is : v = PkPk-1...P0AP0P1...Pkek*/
            for (i = 0; i <= k; i++)
            {
                applyHouseholder(GRID, A, v, Mptr(H, 0, i, mp), i, v);
            }
            /* generate and apply the last transformation */
            generateHouseholder(GRID, A, v, u, k + 1, &tmp);

            /* apply this transformation: v<-Pk+1v */
            HPL_indxg2lp(&index, &pindex, k + 1, A->nb, A->nb, 0, GRID->nprow);
            if (GRID->myrow == pindex)
            {
                v[index] = tmp;
            }
            for (i = 0; i < mp; ++i)
            {
                index = HPL_indxl2g(i, A->nb, A->nb, GRID->myrow, 0, GRID->nprow);
                if (index > k + 1)
                {
                    v[i] = 0;
                }
            }

            /* now v is the kth column of Hm. Apply Givens rotations on v */

            tmp = 0;
            /* generate and apply Givens rotations on v and w */
            /* sinus and cosus store the previous rotation parameters */
            givens_rotations(GRID, A, v, w, R, sinus, cosus, k, MM);
            /* tmp stored the last element of w, which is the current residual */
            tmp = w[k + 1];
            /* store the current error */
            currenterror = fabs(tmp);
            // HPL_barrier(GRID->all_comm);
            // if (GRID->iam == 0)
            // printf("Err: %.16f, start = %d\n", currenterror, start);fflush(stdout);
            // HPL_barrier(GRID->all_comm);
            /* check if the solution is good enough */
            if (currenterror < TOL)
            {
                ready = 1;
                break;
            }
        }
        /* ------------------------------------ */
        /*  Solve for new solution x            */
        /* ------------------------------------ */
        if (k == MM)
        {
            --k;
        }

        // if (GRID->iam == 0)
        // {
        //     printf("Before:\n");
        //     printf("w = :\n");
        //     print_vector(w, MM, 1);
        //     printf("R = :\n");
        //     print_matrix(R, MM, MM, MM, 1);
        // }
        /* solve Ry = w, R is upper-tri, and w will be overwritten by solution y */
        blas::trsv<double, double>(blas::Layout::ColMajor, blas::Uplo::Upper, blas::Op::NoTrans, blas::Diag::NonUnit, k + 1, R, MM, w, 1);
        // if (GRID->iam == 0)
        // {
        //     printf("After:\n");
        //     printf("w = :\n");
        //     print_vector(w, MM, 1);
        // }
        /* calculate the new solution */
        for (i = 0; i <= k; ++i)
        {
            /* load v with unit vector, with v[i] = 1 */
            memset(v, 0, mp * sizeof(double));
            HPL_indxg2lp(&index, &pindex, i, A->nb, A->nb, 0, GRID->nprow);
            if (GRID->myrow == pindex)
            {
                v[index] = 1;
            }
            /* apply last i + 1 householder transformations in reverse order */
            for (j = i; j >= 0; --j)
            {
                applyHouseholder(GRID, A, v, Mptr(H, 0, j, mp), j, v);
            }

            redB2X(GRID, A, v, xt);

            /* update x: perform x += yi*vi */
            for (j = 0; j < nq; ++j)
            {
                x[j] += xt[j] * w[i];
            }
        }

        // /* there is initial guess stored in x here from last iteration */
        // /* calculate v = Ax */
        // blas::gemv<double, double, double>( blas::Layout::ColMajor, blas::Op::NoTrans, mp, nq, HPL_rone,
        //         A->A, A->ld, x, 1, 0, v, 1 );
        // HPL_all_reduce(v, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);

        // /* preconditioning A */
        // if (prec)
        // {
        //     if (GRID->mycol == tarcol)
        //     {
        //         memcpy(bptr, v, mp*sizeof(double));
        //     }
        //     HPL_pLdtrsv(GRID, factors);
        //     redX2B(GRID, factors, factors->X, v);

        //     if (GRID->mycol == tarcol)
        //     {
        //         memcpy(bptr, v, mp*sizeof(double));
        //     }
        //     HPL_pdtrsv(GRID, factors);
        //     redX2B(GRID, factors, factors->X, v);
        // }

        // norm = 0;
        // /* v = rhs - v = rhs - Ax */
        // for (i = 0; i < mp; ++i)
        // {
        //     v[i] = rhs[i] - v[i];
        //     norm += v[i]*v[i];
        // }

        // HPL_all_reduce(&norm, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm);
        // norm = sqrt(norm);

        // HPL_barrier(GRID->all_comm);
        // if (GRID->iam == 0)
        //     printf("currenterror = %.16f, norm = %.16f\n", currenterror, norm);
        // fflush(stdout);
        // HPL_barrier(GRID->all_comm);

        /* if the error is small enough, stop. 
            otherwise another iteration will be initiated. */
        if (currenterror < TOL)
        {
            ready = 1;
        }

        // printf("Restart!\n");
        // fflush(stdout);
    }

    // printf("Final! From process %d\n", GRID->iam);
    // fflush(stdout);
    /* check if we have done maximum number of starts */
    if (start > MAXIT)
    {
        start = MAXIT;
    }

    if (v)
        free(v);
    if (u)
        free(u);
    if (H)
        free(H);
    if (rhs)
        free(rhs);
    if (cosus)
        free(cosus);
    if (sinus)
        free(sinus);
    if (w)
        free(w);
    if (R)
        free(R);
    if (xt)
        free(xt);

    /* return total number of iterations performed */
    return (start * MM + k + 1);

    /* end of HPL_pgmres() */
}

// https://github.com/schuangs/hpl-ai-with-IR/blob/master/src/pir/HPL_pir.c

/*
 * By Junkang Huang, August, 2020
 * 
 * based on the paper:
 *    A New  Analysis  Of Iterative Refinement And Its  Application To 
 *    Accurate Solution Of Ill Conditioned Sparse Linear Systems
 *    ---by Carson, Erin & Higham, Nicholas J., 2017
 */

static void HPL_pir(
    HPL_T_grid *GRID,
    HPLAI_T_palg *ALGO,
    HPL_T_pmat *A,
    HPL_T_pmat *factors,
    double PRE, /* solution tolerance */
    int IR,
    int MM,    /* restart size for GMRES */
    int MAXIT, /* maximum number of GMRES iteration */
    double TOL)
{
    /* 
 * Purpose
 * =======
 *
 * HPL_pir performs iterative refinement procesure to enhance the accur-
 * acy of  the solution  of linear system obtained  by LU factorization. 
 * Parallel  GMRES algorithm  is used  as the inner solver to solve  the 
 * inner correct equation Ad = r.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPL_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * ALGO    (global input)                HPL_T_palg *
 *         On entry,  ALGO  points to  the data structure containing the
 *         algorithmic parameters to be used for this test.
 * 
 * A       (local input/output)          HPL_T_pmat *
 *         On entry, A points to the data structure containing the local
 *         array information. 
 * FA      (local input/output)          HPL_T_pmat *
 *         On entry, A points to the data structure containing the local
 *         array information. It serves as lower precision version if A.
 *
 * ---------------------------------------------------------------------
 */

    /*
 * .. Local Variables ..
 */
    int i, j;
    int mp, nq, n, nb, npcol, nprow, myrow, mycol, tarcol;
    double *Bptr, *res, *d;
    double norm;

    /* ..
 * .. Executable Statements ..
 */
    mp = A->mp;
    nq = A->nq - 1;
    n = A->n;
    nb = A->nb;
    npcol = GRID->npcol;
    nprow = GRID->nprow;
    mycol = GRID->mycol;
    myrow = GRID->myrow;

    /* point to local rhs area */
    Bptr = Mptr(A->A, 0, nq, A->ld);

    /*
 * allocate  space  for  factors, which  contains LU factors of  double
 * precision, which is in fact of lower precision, just stored in double. 
 */
    /*
 * Convert initial solution ( which is obtained through lower precision 
 * LU factorization) into higher precision
 */
    /*
    for (i = 0; i < nq; ++i)
    {
        *(A->X + i) = (double)(*(factors->X + i));
        *(factors->X + i) = 0;
    }
    */
    blas::copy<double, double>(nq, factors->X, 1, A->X, 1);
    //待检查能不能改成上面这个

    /*
 * allocate space for residual vector, correction vectors.
 */
    res = (double *)malloc(mp * sizeof(double));
    d = (double *)malloc(nq * sizeof(double));

    /*
 * tarcol is the process column containing b
 */
    tarcol = HPL_indxg2p(n, nb, nb, 0, npcol);

    /*
 * Iterative Refinement
 */
    for (i = 0; i < IR; ++i)
    {
        /*
        if (GRID->iam == 0)
            printf("IR Loop %d\n", i);
        */
        /* Calculate residual in double precision */
        if (mycol == tarcol)
        {
            memcpy(res, Bptr, mp * sizeof(double));
            blas::gemv<double, double, double>(blas::Layout::ColMajor, blas::Op::NoTrans, mp, nq, -HPL_rone,
                                               A->A, A->ld, A->X, 1, HPL_rone, res, 1);
        }
        else if (nq > 0)
        {
            memset(res, 0, mp * sizeof(double));
            blas::gemv<double, double, double>(blas::Layout::ColMajor, blas::Op::NoTrans, mp, nq, -HPL_rone,
                                               A->A, A->ld, A->X, 1, HPL_rzero, res, 1);
        }
        else
        {
            for (j = 0; j < mp; ++j)
                res[j] = HPL_rzero;
        }

        if (mp > 0)
            HPL_all_reduce(res, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);

        norm = 0;
        for (j = 0; j < mp; ++j)
        {
            norm += res[j] * res[j];
        }
        HPL_all_reduce(&norm, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm);
        norm = sqrt(norm);

        if (norm < PRE)
            break;

        /* 
    * Solve correction  equation using preconditioned  GMRES  method in mix
    * precision.  
    */
        memset(d, 0, nq * sizeof(double));
        HPL_pgmres(GRID, A, factors, res, d, TOL, MM, MAXIT);
        /* 
    * update X with d
    */
        blas::axpy<double, double>(nq, 1, d, 1, A->X, 1);
    }

    /* free dynamic memories */
    if (d)
        free(d);
    if (res)
        free(res);

    /*
 * End of HPL_pir
 */
}

template <typename T1, typename T2>
static void HPLAI_pmat_cpy(
    T1 *DST,
    const T2 *SRC)
{
    DST->n = SRC->n;
    DST->nb = SRC->nb;
    DST->ld = SRC->ld;
    DST->mp = SRC->mp;
    DST->nq = SRC->nq;
    DST->info = SRC->info;
    blas::copy(SRC->nq * SRC->ld, SRC->A, 1, DST->A, 1);
    blas::copy(SRC->nq, SRC->X, 1, DST->X, 1);
}

template <typename T1, typename T2, typename T3>
static void HPLAI_pmat_new(
    T1 *DST,
    const T2 *SRC,
    HPLAI_T_palg *ALGO,
    void **vptr,
    T3 *DSTA)
{
    *vptr = (void *)malloc(
        ((size_t)(ALGO->align) + (size_t)(SRC->ld + 1) * (size_t)(SRC->nq)) * sizeof(DST->A[0]));
    if (*vptr == NULL)
        HPLAI_pabort(__LINE__, "HPLAI_pmat_new", "Memory allocation failed");

#ifdef HPL_CALL_VSIPL
    DST->block = vsip_blockbind_d((vsip_scalar_d *)(SRC->A),
                                  (vsip_length)(SRC->ld * SRC->nq),
                                  VSIP_MEM_NONE);
#endif

    DST->A = (T3 *)HPL_PTR((*vptr), ((size_t)(ALGO->align) * sizeof(DST->A[0])));

    DST->X = Mptr(DST->A, 0, SRC->nq, SRC->ld);

    HPLAI_pmat_cpy(DST, SRC);
}

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    void HPLAI_pdgesv(
        HPL_T_grid *GRID,
        HPLAI_T_palg *ALGO,
        HPL_T_pmat *A)
#else
void HPLAI_pdgesv(GRID, ALGO, A)
    HPL_T_grid *GRID;
HPL_T_palg *ALGO;
HPLAI_T_pmat *A;
#endif
    {
        void *vptr_FA, *vptr_factors;
        HPLAI_T_pmat FA;
        HPL_T_pmat factors;
        HPLAI_pmat_new(&FA, A, ALGO, &vptr_FA, FA.A);

        HPLAI_pagesv(GRID, ALGO, &FA);

#ifdef HPLAI_PMAT_REGEN
        HPLAI_pmat_cpy(A, &FA);
        if (vptr_FA)
            free(vptr_FA);
        HPLAI_pmat_new(&factors, A, ALGO, &vptr_factors, factors.A);
        HPLAI_pdmatgen(GRID, A->n, A->n + 1, A->nb, A->A, A->ld, HPL_ISEED);
#else
    HPLAI_pmat_new(&factors, &FA, ALGO, &vptr_factors, factors.A);
    if (vptr_FA)
        free(vptr_FA);
#endif

        HPL_pir(GRID, ALGO, A, &factors, 1e-14, 1, 50, 1, DBL_EPSILON / 2.0 / ((double)A->n / 4.0));

        if (vptr_factors)
            free(vptr_factors);
    }

#ifdef __cplusplus
}
#endif
