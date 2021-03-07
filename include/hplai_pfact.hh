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
#ifndef HPLAI_PFACT_HH
#define HPLAI_PFACT_HH

#ifdef __cplusplus
extern "C"
{
#endif
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_pfact.h"

#ifdef __cplusplus
}
#endif

#include "hplai_misc.hh"

#include "hplai_panel.hh"
/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */

#ifdef __cplusplus
extern "C"
{
#endif

#define HPLAI_LEFT_LOOKING HPL_LEFT_LOOKING
#define HPLAI_CROUT HPL_CROUT
#define HPLAI_RIGHT_LOOKING HPL_RIGHT_LOOKING
#define HPLAI_T_FACT HPL_T_FACT

#define HPLAI_1RING HPL_1RING
#define HPLAI_1RING_M HPL_1RING_M
#define HPLAI_2RING HPL_2RING
#define HPLAI_2RING_M HPL_2RING_M
#define HPLAI_BLONG HPL_BLONG
#define HPLAI_BLONG_M HPL_BLONG_M
#define HPLAI_T_TOP HPL_T_TOP

typedef struct HPLAI_S_pmat
{
#ifdef HPL_CALL_VSIPL
    vsip_block_d *block;
#endif
    HPLAI_T_AFLOAT *A; /* pointer to local piece of A */
    HPLAI_T_AFLOAT *X; /* pointer to solution vector */
    int n;             /* global problem size */
    int nb;            /* blocking factor */
    int ld;            /* local leading dimension */
    int mp;            /* local number of rows */
    int nq;            /* local number of columns */
    int info;          /* computational flag */
} HPLAI_T_pmat;

struct HPLAI_S_palg;

typedef struct HPLAI_S_panel
{
    struct HPL_S_grid *grid;   /* ptr to the process grid */
    struct HPLAI_S_palg *algo; /* ptr to the algo parameters */
    HPLAI_T_pmat *pmat;        /* ptr to the local array info */
    HPLAI_T_AFLOAT *A;         /* ptr to trailing part of A */
    HPLAI_T_AFLOAT *WORK;      /* work space */
    HPLAI_T_AFLOAT *L2;        /* ptr to L */
    HPLAI_T_AFLOAT *L1;        /* ptr to jb x jb upper block of A */
    HPLAI_T_AFLOAT *DPIV;      /* ptr to replicated jb pivot array */
    HPLAI_T_AFLOAT *DINFO;     /* ptr to replicated scalar info */
    HPLAI_T_AFLOAT *U;         /* ptr to U */
    int *IWORK;                /* integer workspace for swapping */
    void ***buffers[2];        /* buffers for panel bcast */
    int counts[2];             /* counts for panel bcast */
    MPI_Datatype dtypes[2];    /* data types for panel bcast */
    MPI_Request request[1];    /* requests for panel bcast */
    MPI_Status status[1];      /* status for panel bcast */
    int nb;                    /* distribution blocking factor */
    int jb;                    /* panel width */
    int m;                     /* global # of rows of trailing part of A */
    int n;                     /* global # of cols of trailing part of A */
    int ia;                    /* global row index of trailing part of A */
    int ja;                    /* global col index of trailing part of A */
    int mp;                    /* local # of rows of trailing part of A */
    int nq;                    /* local # of cols of trailing part of A */
    int ii;                    /* local row index of trailing part of A */
    int jj;                    /* local col index of trailing part of A */
    int lda;                   /* local leading dim of array A */
    int prow;                  /* proc. row owning 1st row of trail. A */
    int pcol;                  /* proc. col owning 1st col of trail. A */
    int msgid;                 /* message id for panel bcast */
    int ldl2;                  /* local leading dim of array L2 */
    int len;                   /* length of the buffer to broadcast */
#ifdef HPL_CALL_VSIPL
    vsip_block_d *Ablock;  /* A block */
    vsip_block_d *L1block; /* L1 block */
    vsip_block_d *L2block; /* L2 block */
    vsip_block_d *Ublock;  /* U block */
#endif
} HPLAI_T_panel;
/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */
typedef void (*HPLAI_T_PFA_FUN)(HPLAI_T_panel *, const int, const int, const int,
                                HPLAI_T_AFLOAT *);
typedef void (*HPLAI_T_RFA_FUN)(HPLAI_T_panel *, const int, const int, const int,
                                HPLAI_T_AFLOAT *);
typedef void (*HPLAI_T_UPD_FUN)(HPLAI_T_panel *, int *, HPLAI_T_panel *, const int);

#define HPLAI_SWAP00 HPL_SWAP00
#define HPLAI_SWAP01 HPL_SWAP01
#define HPLAI_SW_MIX HPL_SW_MIX
#define HPLAI_NO_SWP HPL_NO_SWP
#define HPLAI_T_SWAP HPL_T_SWAP

typedef struct HPLAI_S_palg
{
    HPLAI_T_TOP btopo;     /* row broadcast topology */
    int depth;             /* look-ahead depth */
    int nbdiv;             /* recursive division factor */
    int nbmin;             /* recursion stopping criterium */
    HPLAI_T_FACT pfact;    /* panel fact variant */
    HPLAI_T_FACT rfact;    /* recursive fact variant */
    HPLAI_T_PFA_FUN pffun; /* panel fact function ptr */
    HPLAI_T_RFA_FUN rffun; /* recursive fact function ptr */
    HPLAI_T_UPD_FUN upfun; /* update function */
    HPLAI_T_SWAP fswap;    /* Swapping algorithm */
    int fsthr;             /* Swapping threshold */
    int equil;             /* Equilibration */
    int align;             /* data alignment constant */
} HPLAI_T_palg;
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void HPLAI_alocmax
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));

void HPLAI_alocswpN
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_alocswpT
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_pamxswp
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));

void HPLAI_papancrN
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_papancrT
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_papanllN
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_papanllT
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_papanrlN
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_papanrlT
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));

void HPLAI_parpancrN
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_parpancrT
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_parpanllN
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_parpanllT
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_parpanrlN
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));
void HPLAI_parpanrlT
    STDC_ARGS((
        HPLAI_T_panel *,
        const int,
        const int,
        const int,
        HPLAI_T_AFLOAT *));

void HPLAI_pafact
    STDC_ARGS((
        HPLAI_T_panel *));

#ifdef __cplusplus
}
#endif

#endif
/*
 * End of hplai_pfact.hh
 */
