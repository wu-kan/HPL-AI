#ifndef HPL_AI_PANEL_H
#define HPL_AI_PANEL_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_ai.h"
#include "hpl_ai_pgesv.h"
/*
 * ---------------------------------------------------------------------
 * Data Structures
 * ---------------------------------------------------------------------
 */
typedef struct HPL_AI_S_panel
{
   struct HPL_S_grid   * grid;             /* ptr to the process grid */
   struct HPL_S_palg   * algo;          /* ptr to the algo parameters */
   HPL_AI_T_pmat       * pmat;         /* ptr to the local array info */
   HPL_AI_T_afloat     * A;              /* ptr to trailing part of A */
   HPL_AI_T_afloat     * WORK;                          /* work space */
   HPL_AI_T_afloat     * L2;                              /* ptr to L */
   HPL_AI_T_afloat     * L1;       /* ptr to jb x jb upper block of A */
   HPL_AI_T_afloat     * DPIV;    /* ptr to replicated jb pivot array */
   HPL_AI_T_afloat     * DINFO;      /* ptr to replicated scalar info */
   HPL_AI_T_afloat     * U;                               /* ptr to U */
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
} HPL_AI_T_panel;

/*
 * ---------------------------------------------------------------------
 * panel function prototypes
 * ---------------------------------------------------------------------
 */

void                             HPL_AI_papanel_new
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   const int,
   const int,
   const int,
   HPL_AI_T_pmat *,
   const int,
   const int,
   const int,
   HPL_AI_T_panel * *
) );
void                             HPL_AI_papanel_init
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   const int,
   const int,
   const int,
   HPL_AI_T_pmat *,
   const int,
   const int,
   const int,
   HPL_AI_T_panel *
) );
int                              HPL_AI_papanel_disp
STDC_ARGS( (
   HPL_AI_T_panel * *
) );
int                              HPL_AI_papanel_free
STDC_ARGS( (
   HPL_AI_T_panel *
) );

#endif
/*
 * End of hpl_ai_panel.h
 */
