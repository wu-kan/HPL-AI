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
#ifndef HPLAI_GRID_HH
#define HPLAI_GRID_HH
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */

#ifdef __cplusplus
extern "C"
{
#endif

#include "hpl_grid.h"

#ifdef __cplusplus
}
#endif

#include "hplai_misc.hh"

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */
#define HPLAI_ROW_MAJOR HPL_ROW_MAJOR
#define HPLAI_COLUMN_MAJOR HPL_COLUMN_MAJOR
#define HPLAI_T_ORDER HPL_T_ORDER
#define HPLAI_S_grid HPL_S_grid
#define HPLAI_T_grid HPL_T_grid

/*
 * ---------------------------------------------------------------------
 * Data Structures
 * ---------------------------------------------------------------------
 */
typedef void (*HPLAI_T_OP_AFLOAT)(const int, const void *, void *);

#define HPLAI_grid_init HPL_grid_init
#define HPLAI_grid_exit HPL_grid_exit
#define HPLAI_grid_info HPL_grid_info
#define HPLAI_pnum HPL_pnum
#define HPLAI_barrier HPL_barrier
int HPLAI_broadcast_AFLOAT
    STDC_ARGS((
        void *,
        const int,
        const int,
        MPI_Comm));
int HPLAI_reduce_AFLOAT
    STDC_ARGS((
        void *,
        const int,
        const HPLAI_T_OP_AFLOAT,
        const int,
        MPI_Comm));
int HPLAI_all_reduce_AFLOAT
    STDC_ARGS((
        void *,
        const int,
        const HPLAI_T_OP_AFLOAT,
        MPI_Comm));

void HPLAI_max_AFLOAT
    STDC_ARGS((
        const int,
        const void *,
        void *));
void HPLAI_min_AFLOAT
    STDC_ARGS((
        const int,
        const void *,
        void *));
void HPLAI_sum_AFLOAT
    STDC_ARGS((
        const int,
        const void *,
        void *));

#ifdef __cplusplus
}
#endif

#endif
/*
 * End of hplai_grid.hh
 */
