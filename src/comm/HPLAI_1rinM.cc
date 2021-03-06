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

#ifdef HPL_NO_MPI_DATATYPE /* The user insists to not use MPI types */
#ifndef HPL_COPY_L         /* and also want to avoid the copy of L ... */
#define HPL_COPY_L         /* well, sorry, can not do that: force the copy */
#endif
#endif

#ifdef STDC_HEADERS
    int HPLAI_binit_1rinM(
        HPLAI_T_panel *PANEL)
#else
int HPLAI_binit_1rinM(PANEL)
    HPLAI_T_panel *PANEL;
#endif
    {
#ifdef HPL_USE_MPI_DATATYPE
        /*
 * .. Local Variables ..
 */
        int ierr;
#endif
        /* ..
 * .. Executable Statements ..
 */
        if (PANEL == NULL)
        {
            return (HPL_SUCCESS);
        }
        if (PANEL->grid->npcol <= 1)
        {
            return (HPL_SUCCESS);
        }
#ifdef HPL_USE_MPI_DATATYPE
#ifdef HPL_COPY_L
        /*
 * Copy the panel into a contiguous buffer
 */
        HPLAI_copyL(PANEL);
#endif
        /*
 * Create the MPI user-defined data type
 */
        ierr = HPLAI_packL(PANEL, 0, PANEL->len, 0);

        return ((ierr == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE));
#else
    /*
 * Force the copy of the panel into a contiguous buffer
 */
    HPLAI_copyL(PANEL);

    return (HPL_SUCCESS);
#endif
    }

#ifdef HPL_USE_MPI_DATATYPE

#define _M_BUFF PANEL->buffers[0]
#define _M_COUNT PANEL->counts[0]
#define _M_TYPE PANEL->dtypes[0]

#else

#define _M_BUFF (void *)(PANEL->L2)
#define _M_COUNT PANEL->len
#define _M_TYPE HPLAI_MPI_AFLOAT

#endif

#ifdef STDC_HEADERS
    int HPLAI_bcast_1rinM(
        HPLAI_T_panel *PANEL,
        int *IFLAG)
#else
int HPLAI_bcast_1rinM(PANEL, IFLAG)
    HPLAI_T_panel *PANEL;
int *IFLAG;
#endif
    {
        /*
 * .. Local Variables ..
 */
        MPI_Comm comm;
        int ierr, go, next, msgid, partner, prev,
            rank, root, size;
        /* ..
 * .. Executable Statements ..
 */
        if (PANEL == NULL)
        {
            *IFLAG = HPL_SUCCESS;
            return (HPL_SUCCESS);
        }
        if ((size = PANEL->grid->npcol) <= 1)
        {
            *IFLAG = HPL_SUCCESS;
            return (HPL_SUCCESS);
        }
        /*
 * Cast phase:  If I am the root process,  then  send message to its two
 * next neighbors. Otherwise, probe for message. If the message is here,
 * then receive it,   and  if I am not the last process of the ring,  or
 * just after the root process, then forward it to the next.  Otherwise,
 * inform the caller that the panel has still not been received.
 */
        rank = PANEL->grid->mycol;
        comm = PANEL->grid->row_comm;
        root = PANEL->pcol;
        msgid = PANEL->msgid;
        next = MModAdd1(rank, size);

        if (rank == root)
        {
            ierr = MPI_Send(_M_BUFF, _M_COUNT, _M_TYPE, next, msgid, comm);
            if ((ierr == MPI_SUCCESS) && (size > 2))
            {
                ierr = MPI_Send(_M_BUFF, _M_COUNT, _M_TYPE, MModAdd1(next, size), msgid, comm);
            }
        }
        else
        {
            prev = MModSub1(rank, size);
            if ((size > 2) &&
                (MModSub1(prev, size) == root))
                partner = root;
            else
                partner = prev;

            ierr = MPI_Iprobe(partner, msgid, comm, &go, &PANEL->status[0]);

            if (ierr == MPI_SUCCESS)
            {
                if (go != 0)
                {
                    ierr = MPI_Recv(_M_BUFF, _M_COUNT, _M_TYPE, partner, msgid,
                                    comm, &PANEL->status[0]);
                    if ((ierr == MPI_SUCCESS) &&
                        (prev != root) && (next != root))
                    {
                        ierr = MPI_Send(_M_BUFF, _M_COUNT, _M_TYPE, next, msgid,
                                        comm);
                    }
                }
                else
                {
                    *IFLAG = HPLAI_KEEP_TESTING;
                    return (*IFLAG);
                }
            }
        }
        /*
 * If the message was received and being forwarded,  return HPL_SUCCESS.
 * If an error occured in an MPI call, return HPL_FAILURE.
 */
        *IFLAG = (ierr == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE);

        return (*IFLAG);
    }

#ifdef STDC_HEADERS
    int HPLAI_bwait_1rinM(
        HPLAI_T_panel *PANEL)
#else
int HPLAI_bwait_1rinM(PANEL)
    HPLAI_T_panel *PANEL;
#endif
    {
#ifdef HPL_USE_MPI_DATATYPE
        /*
 * .. Local Variables ..
 */
        int ierr;
#endif
        /* ..
 * .. Executable Statements ..
 */
        if (PANEL == NULL)
        {
            return (HPL_SUCCESS);
        }
        if (PANEL->grid->npcol <= 1)
        {
            return (HPL_SUCCESS);
        }
/*
 * Release the arrays of request / status / data-types and buffers
 */
#ifdef HPL_USE_MPI_DATATYPE
        ierr = MPI_Type_free(&PANEL->dtypes[0]);
        return ((ierr == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE));
#else
    return (HPL_SUCCESS);
#endif
    }

#ifdef __cplusplus
}
#endif
