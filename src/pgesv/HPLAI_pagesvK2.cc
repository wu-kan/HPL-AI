/*
 * Include files
 */
#include "hplai.hh"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPLAI_pagesvK2
(
   HPLAI_T_grid *                     GRID,
   HPLAI_T_palg *                   ALGO,
   HPLAI_T_pmat *                   A
)
#else
void HPLAI_pagesvK2
( GRID, ALGO, A )
   HPLAI_T_grid *                     GRID;
   HPLAI_T_palg *                   ALGO;
   HPLAI_T_pmat *                   A;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_pagesvK2 factors a N+1-by-N matrix using LU factorization with row
 * partial pivoting.  The main algorithm  is the "right looking" variant
 * with look-ahead.  The  lower  triangular factor is left unpivoted and
 * the pivots are not returned. The right hand side is the N+1 column of
 * the coefficient matrix.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPLAI_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * ALGO    (global input)                HPLAI_T_palg *
 *         On entry,  ALGO  points to  the data structure containing the
 *         algorithmic parameters.
 *
 * A       (local input/output)          HPLAI_T_pmat *
 *         On entry, A points to the data structure containing the local
 *         array information.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   HPLAI_T_panel              * p, * * panel = NULL;
   HPLAI_T_UPD_FUN            HPLAI_paupdate; 
   int                        N, depth, icurcol=0, j, jb, jj=0, jstart,
                              k, mycol, n, nb, nn, npcol, nq,
                              tag=MSGID_BEGIN_FACT, test=HPL_KEEP_TESTING;
#ifdef HPL_PROGRESS_REPORT
   double start_time, time, gflops;
#endif
/* ..
 * .. Executable Statements ..
 */
   mycol = GRID->mycol; npcol        = GRID->npcol;
   depth = ALGO->depth; HPLAI_paupdate = ALGO->upfun;
   N     = A->n;        nb           = A->nb;

   if( N <= 0 ) return;

#ifdef HPL_PROGRESS_REPORT
   start_time = HPL_timer_walltime();
#endif

/*
 * Allocate a panel list of length depth + 1 (depth >= 1)
 */
   panel = (HPLAI_T_panel **)malloc( (size_t)(depth+1) * sizeof( HPLAI_T_panel *) );
   if( panel == NULL )
   { HPL_pabort( __LINE__, "HPLAI_pagesvK2", "Memory allocation failed" ); }
/*
 * Create and initialize the first depth panels
 */
   nq = HPL_numroc( N+1, nb, nb, mycol, 0, npcol ); nn = N; jstart = 0;

   for( k = 0; k < depth; k++ )
   {
      jb = Mmin( nn, nb );
      HPLAI_papanel_new( GRID, ALGO, nn, nn+1, jb, A, jstart, jstart,
                       tag, &panel[k] );
      nn -= jb; jstart += jb;
      if( mycol == icurcol ) { jj += jb; nq -= jb; }
      icurcol = MModAdd1( icurcol, npcol );
      tag     = MNxtMgid( tag, MSGID_BEGIN_FACT, MSGID_END_FACT );
   }
/*
 * Create last depth+1 panel
 */
   HPLAI_papanel_new( GRID, ALGO, nn, nn+1, Mmin( nn, nb ), A, jstart,
                    jstart, tag, &panel[depth] );
   tag = MNxtMgid( tag, MSGID_BEGIN_FACT, MSGID_END_FACT );
/*
 * Initialize the lookahead - Factor jstart columns: panel[0..depth-1]
 */
   for( k = 0, j = 0; k < depth; k++ )
   {
      jb = jstart - j; jb = Mmin( jb, nb ); j += jb;
/*
 * Factor and broadcast k-th panel
 */
      HPLAI_pafact(         panel[k] );
      (void) HPLAI_binit(   panel[k] );
      do
      { (void) HPLAI_bcast( panel[k], &test ); }
      while( test != HPL_SUCCESS );
      (void) HPLAI_bwait(   panel[k] );
/*
 * Partial update of the depth-k-1 panels in front of me
 */
      if( k < depth - 1 )
      {
         nn = HPL_numrocI( jstart-j, j, nb, nb, mycol, 0, npcol );
         HPLAI_paupdate( NULL, NULL, panel[k], nn );
      }
   }
/*
 * Main loop over the remaining columns of A
 */
   for( j = jstart; j < N; j += nb )
   {
      n = N - j; jb = Mmin( n, nb );
#ifdef HPL_PROGRESS_REPORT
      /* if this is process 0,0 and not the first panel */
      if ( GRID->myrow == 0 && mycol == 0 && j > 0 ) 
      {
          time = HPL_timer_walltime() - start_time;
          gflops = 2.0*(N*(double)N*N - n*(double)n*n)/3.0/(time > 0.0 ? time : 1e-6)/1e9;
          HPL_fprintf( stdout, "Column=%09d Fraction=%4.1f%% Gflops=%9.3e\n", j, j*100.0/N, gflops);
      }
#endif
/*
 * Initialize current panel - Finish latest update, Factor and broadcast
 * current panel
 */
      (void) HPLAI_papanel_free( panel[depth] );
      HPLAI_papanel_init( GRID, ALGO, n, n+1, jb, A, j, j, tag, panel[depth] );

      if( mycol == icurcol )
      {
         nn = HPL_numrocI( jb, j, nb, nb, mycol, 0, npcol );
         for( k = 0; k < depth; k++ )   /* partial updates 0..depth-1 */
            (void) HPLAI_paupdate( NULL, NULL, panel[k], nn );
         HPLAI_pafact(       panel[depth] );    /* factor current panel */
      }
      else { nn = 0; }
          /* Finish the latest update and broadcast the current panel */
      (void) HPLAI_binit( panel[depth] );
      HPLAI_paupdate( panel[depth], &test, panel[0], nq-nn );
      (void) HPLAI_bwait( panel[depth] );
/*
 * Circular  of the panel pointers:
 * xtmp = x[0]; for( k=0; k < depth; k++ ) x[k] = x[k+1]; x[d] = xtmp;
 *
 * Go to next process row and column - update the message ids for broadcast
 */
      p = panel[0]; for( k = 0; k < depth; k++ ) panel[k] = panel[k+1];
      panel[depth] = p;

      if( mycol == icurcol ) { jj += jb; nq -= jb; }
      icurcol = MModAdd1( icurcol, npcol );
      tag     = MNxtMgid( tag, MSGID_BEGIN_FACT, MSGID_END_FACT );
   }
/*
 * Clean-up: Finish updates - release panels and panel list
 */
   nn = HPL_numrocI( 1, N, nb, nb, mycol, 0, npcol );
   for( k = 0; k < depth; k++ )
   {
      (void) HPLAI_paupdate( NULL, NULL, panel[k], nn );
      (void) HPLAI_papanel_disp(  &panel[k] );
   }
   (void) HPLAI_papanel_disp( &panel[depth] );

   if( panel ) free( panel );
/*
 * End of HPLAI_pagesvK2
 */
}

#ifdef __cplusplus
}
#endif
