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
    int HPLAI_reduce_AFLOAT(
        void *BUFFER,
        const int COUNT,
        const HPLAI_T_OP_AFLOAT OP,
        const int ROOT,
        MPI_Comm COMM)
#else
int HPLAI_reduce_AFLOAT(BUFFER, COUNT, OP, ROOT, COMM) void *BUFFER;
const int COUNT;
const HPLAI_T_OP_AFLOAT OP;
const int ROOT;
MPI_Comm COMM;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_reduce_AFLOAT performs a global reduce operation across all processes of
 * a group.  Note that the input buffer is  used as workarray and in all
 * processes but the accumulating process corrupting the original data.
 *
 * Arguments
 * =========
 *
 * BUFFER  (local input/output)          void *
 *         On entry,  BUFFER  points to  the  buffer to be  reduced.  On
 *         exit,  and  in process of rank  ROOT  this array contains the
 *         reduced data.  This  buffer  is also used as workspace during
 *         the operation in the other processes of the group.
 *
 * COUNT   (global input)                const int
 *         On entry,  COUNT  indicates the number of entries in  BUFFER.
 *         COUNT must be at least zero.
 *
 * OP      (global input)                const HPL_T_OP 
 *         On entry, OP is a pointer to the local combine function.
 *
 * ROOT    (global input)                const int
 *         On entry, ROOT is the coordinate of the accumulating process.
 *
 * COMM    (global/local input)          MPI_Comm
 *         The MPI communicator identifying the process collection.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        MPI_Status status;
        void *buffer = NULL;
        int hplerr = MPI_SUCCESS, d = 1, i, ip2 = 1, mask = 0,
            mpierr, mydist, partner, rank, size,
            tag = MSGID_BEGIN_COLL;
        /* ..
 * .. Executable Statements ..
 */
        if (COUNT <= 0)
            return (MPI_SUCCESS);
        mpierr = MPI_Comm_size(COMM, &size);
        if (size == 1)
            return (MPI_SUCCESS);
        mpierr = MPI_Comm_rank(COMM, &rank);
        i = size - 1;
        while (i > 1)
        {
            i >>= 1;
            d++;
        }

        buffer = (void *)((HPLAI_T_AFLOAT *)malloc((size_t)(COUNT) *
                                                   sizeof(HPLAI_T_AFLOAT)));

        if (!(buffer))
        {
            HPL_pabort(__LINE__, "HPLAI_reduce_AFLOAT", "Memory allocation failed");
        }

        if ((mydist = MModSub(rank, ROOT, size)) == 0)
        {
            do
            {
                mpierr = MPI_Recv(buffer, COUNT, HPLAI_MPI_AFLOAT,
                                  MModAdd(ROOT, ip2, size), tag, COMM,
                                  &status);
                if (mpierr != MPI_SUCCESS)
                    hplerr = mpierr;
                OP(COUNT, buffer, BUFFER);
                ip2 <<= 1;
                d--;
            } while (d);
        }
        else
        {
            do
            {
                if ((mydist & mask) == 0)
                {
                    partner = mydist ^ ip2;

                    if (mydist & ip2)
                    {
                        partner = MModAdd(ROOT, partner, size);
                        mpierr = MPI_Send(BUFFER, COUNT, HPLAI_MPI_AFLOAT,
                                          partner, tag, COMM);
                    }
                    else if (partner < size)
                    {
                        partner = MModAdd(ROOT, partner, size);
                        mpierr = MPI_Recv(buffer, COUNT, HPLAI_MPI_AFLOAT,
                                          partner, tag, COMM, &status);
                        OP(COUNT, buffer, BUFFER);
                    }
                    if (mpierr != MPI_SUCCESS)
                        hplerr = mpierr;
                }
                mask ^= ip2;
                ip2 <<= 1;
                d--;
            } while (d);
        }
        if (buffer)
            free(buffer);

        return (hplerr);
        /*
 * End of HPLAI_reduce_AFLOAT
 */
    }

#ifdef __cplusplus
}
#endif