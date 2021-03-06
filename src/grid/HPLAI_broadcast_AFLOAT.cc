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
    int HPLAI_broadcast_AFLOAT(
        void *BUFFER,
        const int COUNT,
        const int ROOT,
        MPI_Comm COMM)
#else
int HPLAI_broadcast_AFLOAT(BUFFER, COUNT, ROOT, COMM) void *BUFFER;
const int COUNT;
const int ROOT;
MPI_Comm COMM;
#endif
    {
        if (COUNT > 0)
            return MPI_Bcast(BUFFER, COUNT, HPLAI_MPI_AFLOAT, ROOT, COMM);
        /*
 * .. Local Variables ..
 */
        int hplerr = MPI_SUCCESS, ip2 = 1, kk, mask = 1,
            mpierr, mydist, partner, rank, size,
            tag = MSGID_BEGIN_COLL;
        MPI_Status status;
        /* ..
 * .. Executable Statements ..
 */
        if (COUNT <= 0)
            return (MPI_SUCCESS);
        mpierr = MPI_Comm_size(COMM, &size);
        if (size <= 1)
            return (mpierr);
        mpierr = MPI_Comm_rank(COMM, &rank);

        kk = size - 1;
        while (kk > 1)
        {
            kk >>= 1;
            ip2 <<= 1;
            mask <<= 1;
            mask++;
        }
        mydist = MModSub(rank, ROOT, size);

        do
        {
            mask ^= ip2;
            if ((mydist & mask) == 0)
            {
                partner = mydist ^ ip2;

                if (mydist & ip2)
                {
                    partner = MModAdd(ROOT, partner, size);
                    mpierr = MPI_Recv(BUFFER, COUNT, HPLAI_MPI_AFLOAT,
                                      partner, tag, COMM, &status);
                }
                else if (partner < size)
                {
                    partner = MModAdd(ROOT, partner, size);
                    mpierr = MPI_Send(BUFFER, COUNT, HPLAI_MPI_AFLOAT,
                                      partner, tag, COMM);
                }
                if (mpierr != MPI_SUCCESS)
                    hplerr = mpierr;
            }
            ip2 >>= 1;
        } while (ip2);

        return (hplerr);
    }

#ifdef __cplusplus
}
#endif
