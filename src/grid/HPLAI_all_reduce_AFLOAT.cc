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
    int HPLAI_all_reduce_AFLOAT(
        void *BUFFER,
        const int COUNT,
        const HPLAI_T_OP_AFLOAT OP,
        MPI_Comm COMM)
#else
int HPLAI_all_reduce_AFLOAT(BUFFER, COUNT, DTYPE, OP, COMM) void *BUFFER;
const int COUNT;
const HPLAI_T_OP_AFLOAT OP;
MPI_Comm COMM;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * HPLAI_all_reduce_AFLOAT performs   a   global   reduce  operation  across  all
 * processes of a group leaving the results on all processes.
 *
 * Arguments
 * =========
 *
 * BUFFER  (local input/global output)   void *
 *         On entry,  BUFFER  points to  the  buffer to be combined.  On
 *         exit, this array contains the combined data and  is identical
 *         on all processes in the group.
 *
 * COUNT   (global input)                const int
 *         On entry,  COUNT  indicates the number of entries in  BUFFER.
 *         COUNT must be at least zero.
 *
 * OP      (global input)                const HPLAI_T_OP_AFLOAT 
 *         On entry, OP is a pointer to the local combine function.
 *
 * COMM    (global/local input)          MPI_Comm
 *         The MPI communicator identifying the process collection.
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        int hplerr;
        /* ..
 * .. Executable Statements ..
 */
        hplerr = HPLAI_reduce_AFLOAT(BUFFER, COUNT, OP, 0, COMM);
        if (hplerr != MPI_SUCCESS)
            return (hplerr);
        return (HPLAI_broadcast_AFLOAT(BUFFER, COUNT, 0, COMM));
        /*
 * End of HPLAI_all_reduce_AFLOAT
 */
    }

#ifdef __cplusplus
}
#endif
