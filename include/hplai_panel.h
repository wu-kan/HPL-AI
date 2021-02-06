#ifndef HPLAI_PANEL_H
#define HPLAI_PANEL_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_panel.h"
#include "hplai_misc.h"
#include "hplai_grid.h"
/*
 * ---------------------------------------------------------------------
 * Data Structures
 * ---------------------------------------------------------------------
 */
typedef struct HPLAI_S_pmat
{
#ifdef HPL_CALL_VSIPL
   vsip_block_d        * block;
#endif
   HPLAI_T_AFLOAT     * A;            /* pointer to local piece of A */
   HPLAI_T_AFLOAT     * X;             /* pointer to solution vector */
   int                 n;                      /* global problem size */
   int                 nb;                         /* blocking factor */
   int                 ld;                 /* local leading dimension */
   int                 mp;                    /* local number of rows */
   int                 nq;                 /* local number of columns */
   int                 info;                    /* computational flag */
} HPLAI_T_pmat;

typedef struct HPLAI_S_panel
{
   struct HPL_S_grid   * grid;             /* ptr to the process grid */
   struct HPL_S_palg   * algo;          /* ptr to the algo parameters */
   HPLAI_T_pmat       * pmat;         /* ptr to the local array info */
   HPLAI_T_AFLOAT     * A;              /* ptr to trailing part of A */
   HPLAI_T_AFLOAT     * WORK;                          /* work space */
   HPLAI_T_AFLOAT     * L2;                              /* ptr to L */
   HPLAI_T_AFLOAT     * L1;       /* ptr to jb x jb upper block of A */
   HPLAI_T_AFLOAT     * DPIV;    /* ptr to replicated jb pivot array */
   HPLAI_T_AFLOAT     * DINFO;      /* ptr to replicated scalar info */
   HPLAI_T_AFLOAT     * U;                               /* ptr to U */
   int                 * IWORK;     /* integer workspace for swapping */
   void                * * * buffers[2];   /* buffers for panel bcast */
   int                 counts [2];          /* counts for panel bcast */
   MPI_Datatype        dtypes [2];      /* data types for panel bcast */
   MPI_Request         request[1];        /* requests for panel bcast */
   MPI_Status          status [1];          /* status for panel bcast */
   int                 nb;            /* distribution blocking factor */
   int                 jb;                             /* panel width */
   int                 m;   /* global # of rows of trailing part of A */
   int                 n;   /* global # of cols of trailing part of A */
   int                 ia;  /* global row index of trailing part of A */
   int                 ja;  /* global col index of trailing part of A */
   int                 mp;   /* local # of rows of trailing part of A */
   int                 nq;   /* local # of cols of trailing part of A */
   int                 ii;   /* local row index of trailing part of A */
   int                 jj;   /* local col index of trailing part of A */
   int                 lda;           /* local leading dim of array A */
   int                 prow;  /* proc. row owning 1st row of trail. A */
   int                 pcol;  /* proc. col owning 1st col of trail. A */
   int                 msgid;           /* message id for panel bcast */
   int                 ldl2;         /* local leading dim of array L2 */
   int                 len;      /* length of the buffer to broadcast */
#ifdef HPL_CALL_VSIPL
   vsip_block_d        * Ablock;                           /* A block */
   vsip_block_d        * L1block;                         /* L1 block */
   vsip_block_d        * L2block;                         /* L2 block */
   vsip_block_d        * Ublock;                           /* U block */
#endif
} HPLAI_T_panel;

/*
 * ---------------------------------------------------------------------
 * panel function prototypes
 * ---------------------------------------------------------------------
 */

void                             HPLAI_papanel_new
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   const int,
   const int,
   const int,
   HPLAI_T_pmat *,
   const int,
   const int,
   const int,
   HPLAI_T_panel * *
) );
void                             HPLAI_papanel_init
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   const int,
   const int,
   const int,
   HPLAI_T_pmat *,
   const int,
   const int,
   const int,
   HPLAI_T_panel *
) );
int                              HPLAI_papanel_disp
STDC_ARGS( (
   HPLAI_T_panel * *
) );
int                              HPLAI_papanel_free
STDC_ARGS( (
   HPLAI_T_panel *
) );

#endif
/*
 * End of hplai_panel.h
 */
