/*
 * Include files
 */
#include "hplai.hh"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
void HPLAI_patrsv
(
   HPLAI_T_grid *                     GRID,
   HPLAI_T_pmat *                     AMAT
)
#else
void HPLAI_patrsv
( GRID, AMAT )
   HPLAI_T_grid *                     GRID;
   HPLAI_T_pmat *                     AMAT;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPLAI_patrsv solves an upper triangular system of linear equations.
 *  
 * The rhs is the last column of the N by N+1 matrix A. The solve starts
 * in the process  column owning the  Nth  column of A, so the rhs b may
 * need to be moved one process column to the left at the beginning. The
 * routine therefore needs  a column  vector in every process column but
 * the one owning  b. The result is  replicated in all process rows, and
 * returned in XR, i.e. XR is of size nq = LOCq( N ) in all processes.
 *  
 * The algorithm uses decreasing one-ring broadcast in process rows  and
 * columns  implemented  in terms of  synchronous communication point to
 * point primitives.  The  lookahead of depth 1 is used to minimize  the
 * critical path. This entire operation is essentially ``latency'' bound
 * and an estimate of its running time is given by:
 *  
 *    (move rhs) lat + N / ( P bdwth ) +            
 *    (solve)    ((N / NB)-1) 2 (lat + NB / bdwth) +
 *               gam2 N^2 / ( P Q ),                
 *  
 * where  gam2   is an estimate of the   Level 2 BLAS rate of execution.
 * There are  N / NB  diagonal blocks. One must exchange  2  messages of
 * length NB to compute the next  NB  entries of the vector solution, as
 * well as performing a total of N^2 floating point operations.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPLAI_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * AMAT    (local input/output)          HPLAI_T_pmat *
 *         On entry,  AMAT  points  to the data structure containing the
 *         local array information.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   MPI_Comm                   Ccomm, Rcomm;
   HPLAI_T_AFLOAT                     * A=NULL, * Aprev=NULL, * Aptr, * XC=NULL,
                              * XR=NULL, * Xd=NULL, * Xdprev=NULL,
                              * W=NULL;
   int                        Alcol, Alrow, Anpprev, Anp, Anq, Bcol,
                              Cmsgid, GridIsNotPx1, GridIsNot1xQ, Rmsgid,
                              Wfr=0, colprev, kb, kbprev, lda, mycol,
                              myrow, n, n1, n1p, n1pprev=0, nb, npcol,
                              nprow, rowprev, tmp1, tmp2;
/* ..
 * .. Executable Statements ..
 */
#ifdef HPL_DETAILED_TIMING
   HPL_ptimer( HPL_TIMING_PTRSV );
#endif
   if( ( n = AMAT->n ) <= 0 ) return;
   nb = AMAT->nb; lda = AMAT->ld; A = AMAT->A; XR = AMAT->X;

   (void) HPLAI_grid_info( GRID, &nprow, &npcol, &myrow, &mycol );
   Rcomm = GRID->row_comm; Rmsgid = MSGID_BEGIN_PTRSV;
   Ccomm = GRID->col_comm; Cmsgid = MSGID_BEGIN_PTRSV + 1;
   GridIsNot1xQ = ( nprow > 1 ); GridIsNotPx1 = ( npcol > 1 );
/*
 * Move the rhs in the process column owning the last column of A.
 */
   Mnumroc( Anp, n, nb, nb, myrow, 0, nprow );
   Mnumroc( Anq, n, nb, nb, mycol, 0, npcol );

   tmp1  = ( n - 1 ) / nb;
   Alrow = tmp1 - ( tmp1 / nprow ) * nprow;
   Alcol = tmp1 - ( tmp1 / npcol ) * npcol;
   kb    = n    - tmp1 * nb;

   Aptr = (HPLAI_T_AFLOAT *)(A); XC = Mptr( Aptr, 0, Anq, lda );
   Mindxg2p( n, nb, nb, Bcol, 0, npcol );

   if( ( Anp > 0 ) && ( Alcol != Bcol ) )
   {
      if( mycol == Bcol  )
      { (void) HPLAI_send( XC, Anp, Alcol, Rmsgid, Rcomm ); }
      else if( mycol == Alcol )
      { (void) HPLAI_recv( XC, Anp, Bcol,  Rmsgid, Rcomm ); }
   }
   Rmsgid = ( Rmsgid + 2 >
              MSGID_END_PTRSV ? MSGID_BEGIN_PTRSV : Rmsgid + 2 );
   if( mycol != Alcol )
   { for( tmp1=0; tmp1 < Anp; tmp1++ ) XC[tmp1] = HPLAI_rzero; }
/*
 * Set up lookahead
 */
   n1 = ( npcol - 1 ) * nb; n1 = Mmax( n1, nb );
   if( Anp > 0 )
   {
      W = (HPLAI_T_AFLOAT*)malloc( (size_t)(Mmin( n1, Anp )) * sizeof( HPLAI_T_AFLOAT ) );
      if( W == NULL )
      { HPLAI_pabort( __LINE__, "HPLAI_patrsv", "Memory allocation failed" ); }
      Wfr = 1;
   }

   Anpprev = Anp; Xdprev = XR; Aprev = Aptr = Mptr( Aptr, 0, Anq, lda );
   tmp1    = n - kb; tmp1 -= ( tmp2 = Mmin( tmp1, n1 ) );
   MnumrocI( n1pprev, tmp2, Mmax( 0, tmp1 ), nb, nb, myrow, 0, nprow );

   if( myrow == Alrow ) { Anpprev = ( Anp -= kb ); }
   if( mycol == Alcol )
   {
      Aprev = ( Aptr -= lda * kb ); Anq -= kb; Xdprev = ( Xd = XR + Anq );
      if( myrow == Alrow )
      {
         blas::trsv( blas::Layout::ColMajor, blas::Uplo::Upper, blas::Op::NoTrans, blas::Diag::NonUnit,
                    kb, Aptr+Anp, lda, XC+Anp, 1 );
         blas::copy( kb, XC+Anp, 1, Xd, 1 );
      }
   }

   rowprev = Alrow; Alrow = MModSub1( Alrow, nprow );
   colprev = Alcol; Alcol = MModSub1( Alcol, npcol );
   kbprev  = kb; n -= kb;
   tmp1    = n - ( kb = nb ); tmp1 -= ( tmp2 = Mmin( tmp1, n1 ) );
   MnumrocI( n1p, tmp2, Mmax( 0, tmp1 ), nb, nb, myrow, 0, nprow );
/*
 * Start the operations
 */
   while( n > 0 )
   {
      if( mycol == Alcol ) { Aptr -= lda * kb; Anq -= kb; Xd = XR + Anq; }
      if( myrow == Alrow ) { Anp -= kb; }
/*
 * Broadcast  (decreasing-ring)  of  previous solution block in previous
 * process column,  compute  partial update of current block and send it
 * to current process column.
 */
      if( mycol == colprev )
      {
/*
 * Send previous solution block in process row above
 */
         if( myrow == rowprev )
         {
            if( GridIsNot1xQ )
               (void) HPLAI_send( Xdprev, kbprev, MModSub1( myrow, nprow ),
                                Cmsgid, Ccomm );
         }
         else
         {
            (void) HPLAI_recv( Xdprev, kbprev, MModAdd1( myrow, nprow ),
                             Cmsgid, Ccomm );
         } 
/*
 * Compute partial update of previous solution block and send it to cur-
 * rent column
 */
         if( n1pprev > 0 )
         {
            tmp1 = Anpprev - n1pprev;
            blas::gemv( blas::Layout::ColMajor, blas::Op::NoTrans, n1pprev, kbprev,
                       -HPLAI_rone, Aprev+tmp1, lda, Xdprev, 1, HPLAI_rone,
                       XC+tmp1, 1 );
            if( GridIsNotPx1 )
               (void) HPLAI_send( XC+tmp1, n1pprev, Alcol, Rmsgid, Rcomm );
         }
/*
 * Finish  the (decreasing-ring) broadcast of the solution block in pre-
 * vious process column
 */
         if( ( myrow != rowprev ) &&
             ( myrow != MModAdd1( rowprev, nprow ) ) )
            (void) HPLAI_send( Xdprev, kbprev, MModSub1( myrow, nprow ),
                             Cmsgid, Ccomm );
      }
      else if( mycol == Alcol )
      {
/*
 * Current  column  receives  and accumulates partial update of previous
 * solution block
 */
         if( n1pprev > 0 )
         {
            (void) HPLAI_recv( W, n1pprev, colprev, Rmsgid, Rcomm );
            blas::axpy( n1pprev, HPLAI_rone, W, 1, XC+Anpprev-n1pprev, 1 );
         }
      }
/*
 * Solve current diagonal block 
 */
      if( ( mycol == Alcol ) && ( myrow == Alrow ) )
      {
         blas::trsv( blas::Layout::ColMajor, blas::Uplo::Upper, blas::Op::NoTrans, blas::Diag::NonUnit,
                    kb, Aptr+Anp, lda, XC+Anp, 1 );
         blas::copy( kb, XC+Anp, 1, XR+Anq, 1 );
      }
/*
*  Finish previous update
*/
      if( ( mycol == colprev ) && ( ( tmp1 = Anpprev - n1pprev ) > 0 ) )
         blas::gemv( blas::Layout::ColMajor, blas::Op::NoTrans, tmp1, kbprev, -HPLAI_rone,
                    Aprev, lda, Xdprev, 1, HPLAI_rone, XC, 1 );
/*
*  Save info of current step and update info for the next step
*/
      if( mycol == Alcol ) { Xdprev   = Xd; Aprev = Aptr; }
      if( myrow == Alrow ) { Anpprev -= kb; }
      rowprev = Alrow; colprev = Alcol;
      n1pprev = n1p;   kbprev  = kb; n -= kb;
      Alrow = MModSub1( Alrow, nprow ); Alcol = MModSub1( Alcol, npcol );
      tmp1  = n - ( kb = nb ); tmp1 -= ( tmp2 = Mmin( tmp1, n1 ) );
      MnumrocI( n1p, tmp2, Mmax( 0, tmp1 ), nb, nb, myrow, 0, nprow );

      Rmsgid = ( Rmsgid+2 > MSGID_END_PTRSV ? 
                 MSGID_BEGIN_PTRSV   : Rmsgid+2 );
      Cmsgid = ( Cmsgid+2 > MSGID_END_PTRSV ?
                 MSGID_BEGIN_PTRSV+1 : Cmsgid+2 );
   }
/*
 * Replicate last solution block
 */
   if( mycol == colprev )
      (void) HPLAI_broadcast_AFLOAT( (void *)(XR), kbprev, rowprev,
                            Ccomm );

   if( Wfr  ) free( W  );
#ifdef HPL_DETAILED_TIMING
   HPL_ptimer( HPL_TIMING_PTRSV );
#endif
/*
 * End of HPLAI_patrsv
 */
}

#ifdef __cplusplus
}
#endif
