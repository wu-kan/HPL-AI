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

#define I_SEND 0
#define I_RECV 1

#ifdef STDC_HEADERS
void HPLAI_rollN(
    HPLAI_T_panel *PBCST,
    int *IFLAG,
    HPLAI_T_panel *PANEL,
    const int N,
    HPLAI_T_AFLOAT *U,
    const int LDU,
    const int *IPLEN,
    const int *IPMAP,
    const int *IPMAPM1)
#else
void HPLAI_rollN(PBCST, IFLAG, PANEL, N, U, LDU, IPLEN, IPMAP, IPMAPM1)
    HPLAI_T_panel *PBCST;
int *IFLAG;
HPLAI_T_panel *PANEL;
const int N;
HPLAI_T_AFLOAT *U;
const int LDU;
const int *IPLEN;
const int *IPMAP;
const int *IPMAPM1;
#endif
{
    /* 
 * Purpose
 * =======
 *
 * HPLAI_rollN rolls the local arrays containing the local pieces of U, so
 * that on exit to this function  U  is replicated in every process row.
 * In addition, this function probe for the presence of the column panel
 * and forwards it when available.
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
 *         panel (to be rolled) information.
 *
 * N       (local input)                 const int
 *         On entry, N specifies the number of columns of  U.  N must be
 *         at least zero.
 *
 * U       (local input/output)          HPLAI_T_AFLOAT *
 *         On entry,  U  is an array of dimension (LDU,*) containing the
 *         local pieces of U in each process row.
 *
 * LDU     (local input)                 const int
 *         On entry, LDU specifies the local leading dimension of U. LDU
 *         should be at least  MAX(1,IPLEN[NPROW]).
 *
 * IPLEN   (global input)                const int *
 *         On entry, IPLEN is an array of dimension NPROW+1.  This array
 *         is such that IPLEN[i+1] - IPLEN[i] is the number of rows of U
 *         in each process row.
 *
 * IPMAP   (global input)                const int *
 *         On entry, IMAP  is an array of dimension  NPROW.  This  array
 *         contains  the  logarithmic mapping of the processes. In other
 *         words,  IMAP[myrow]  is the absolute coordinate of the sorted
 *         process.
 *
 * IPMAPM1 (global input)                const int *
 *         On entry,  IMAPM1  is an array of dimension NPROW. This array
 *         contains  the inverse of the logarithmic mapping contained in
 *         IMAP: For i in [0.. NPROW) IMAPM1[IMAP[i]] = i.
 *
 * ---------------------------------------------------------------------
 */
    /*
 * .. Local Variables ..
 */
    MPI_Datatype type[2];
    MPI_Status status;
    MPI_Request request;
    MPI_Comm comm;
    int Cmsgid = MSGID_BEGIN_PFACT, ibufR, ibufS,
        ierr = MPI_SUCCESS, il, k, l, lengthR,
        lengthS, mydist, myrow, next, npm1, nprow,
        partner, prev;
    /* ..
 * .. Executable Statements ..
 */
    if (N <= 0)
        return;

    npm1 = (nprow = PANEL->grid->nprow) - 1;
    myrow = PANEL->grid->myrow;
    comm = PANEL->grid->col_comm;
    /*
 * Rolling phase
 */
    mydist = IPMAPM1[myrow];
    prev = IPMAP[MModSub1(mydist, nprow)];
    next = IPMAP[MModAdd1(mydist, nprow)];

    for (k = 0; k < npm1; k++)
    {
        l = (int)((unsigned int)(k) >> 1);

        if (((mydist + k) & 1) != 0)
        {
            il = MModAdd(mydist, l, nprow);
            lengthS = IPLEN[il + 1] - (ibufS = IPLEN[il]);
            il = MModSub(mydist, l + 1, nprow);
            lengthR = IPLEN[il + 1] - (ibufR = IPLEN[il]);
            partner = prev;
        }
        else
        {
            il = MModSub(mydist, l, nprow);
            lengthS = IPLEN[il + 1] - (ibufS = IPLEN[il]);
            il = MModAdd(mydist, l + 1, nprow);
            lengthR = IPLEN[il + 1] - (ibufR = IPLEN[il]);
            partner = next;
        }

        if (lengthR > 0)
        {
            if (ierr == MPI_SUCCESS)
                ierr = MPI_Type_vector(N, lengthR, LDU, HPLAI_MPI_AFLOAT,
                                       &type[I_RECV]);
            if (ierr == MPI_SUCCESS)
                ierr = MPI_Type_commit(&type[I_RECV]);
            if (ierr == MPI_SUCCESS)
                ierr = MPI_Irecv(Mptr(U, ibufR, 0, LDU), 1, type[I_RECV],
                                 partner, Cmsgid, comm, &request);
        }

        if (lengthS > 0)
        {
            if (ierr == MPI_SUCCESS)
                ierr = MPI_Type_vector(N, lengthS, LDU, HPLAI_MPI_AFLOAT,
                                       &type[I_SEND]);
            if (ierr == MPI_SUCCESS)
                ierr = MPI_Type_commit(&type[I_SEND]);
            if (ierr == MPI_SUCCESS)
                ierr = MPI_Send(Mptr(U, ibufS, 0, LDU), 1, type[I_SEND],
                                partner, Cmsgid, comm);
            if (ierr == MPI_SUCCESS)
                ierr = MPI_Type_free(&type[I_SEND]);
        }

        if (lengthR > 0)
        {
            if (ierr == MPI_SUCCESS)
                ierr = MPI_Wait(&request, &status);
            if (ierr == MPI_SUCCESS)
                ierr = MPI_Type_free(&type[I_RECV]);
        }
        /*
 * Probe for column panel - forward it when available
 */
        if (*IFLAG == HPLAI_KEEP_TESTING)
            (void)HPLAI_bcast(PBCST, IFLAG);
    }

    if (ierr != MPI_SUCCESS)
    {
        HPLAI_pabort(__LINE__, "HPLAI_rollN", "MPI call failed");
    }
    /*
 * End of HPLAI_rollN
 */
}
