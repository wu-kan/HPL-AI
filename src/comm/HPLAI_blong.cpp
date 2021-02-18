/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef HPL_NO_MPI_DATATYPE  /* The user insists to not use MPI types */
#ifndef HPL_COPY_L       /* and also want to avoid the copy of L ... */
#define HPL_COPY_L   /* well, sorry, can not do that: force the copy */
#endif
#endif

#define   I_SEND    0
#define   I_RECV    1

#ifdef STDC_HEADERS
int HPLAI_binit_blong
(
   HPLAI_T_panel *              PANEL
)
#else
int HPLAI_binit_blong( PANEL )
   HPLAI_T_panel *              PANEL;
#endif
{
/* ..
 * .. Executable Statements ..
 */
   if( PANEL == NULL )           { return( HPL_SUCCESS ); }
   if( PANEL->grid->npcol <= 1 ) { return( HPL_SUCCESS ); }
#ifdef HPL_USE_MPI_DATATYPE
#ifdef HPL_COPY_L
/*
 * Copy the panel into a contiguous buffer
 */
   HPLAI_copyL( PANEL );
#endif
#else
/*
 * Force the copy of the panel into a contiguous buffer
 */
   HPLAI_copyL( PANEL );
#endif
   return( HPL_SUCCESS );
}

#ifdef HPL_USE_MPI_DATATYPE
 
#define   _M_BUFF_S         PANEL->buffers[I_SEND]
#define   _M_COUNT_S        PANEL->counts[I_SEND]
#define   _M_TYPE_S         PANEL->dtypes[I_SEND]
 
#define   _M_BUFF_R         PANEL->buffers[I_RECV]
#define   _M_COUNT_R        PANEL->counts[I_RECV]
#define   _M_TYPE_R         PANEL->dtypes[I_RECV]
 
#define   _M_ROLL_BUFF_S    PANEL->buffers[I_SEND]
#define   _M_ROLL_COUNT_S   PANEL->counts[I_SEND]
#define   _M_ROLL_TYPE_S    PANEL->dtypes[I_SEND]
 
#define   _M_ROLL_BUFF_R    PANEL->buffers[I_RECV]
#define   _M_ROLL_COUNT_R   PANEL->counts[I_RECV]
#define   _M_ROLL_TYPE_R    PANEL->dtypes[I_RECV]
 
#else
 
#define   _M_BUFF_S         (void *)(PANEL->L2 + ibuf)
#define   _M_COUNT_S        lbuf
#define   _M_TYPE_S         HPLAI_MPI_AFLOAT
 
#define   _M_BUFF_R         (void *)(PANEL->L2 + ibuf)
#define   _M_COUNT_R        lbuf
#define   _M_TYPE_R         HPLAI_MPI_AFLOAT
 
#define   _M_ROLL_BUFF_S    (void *)(PANEL->L2 + ibufS)
#define   _M_ROLL_COUNT_S   lbufS
#define   _M_ROLL_TYPE_S    HPLAI_MPI_AFLOAT
 
#define   _M_ROLL_BUFF_R    (void *)(PANEL->L2 + ibufR)
#define   _M_ROLL_COUNT_R   lbufR
#define   _M_ROLL_TYPE_R    HPLAI_MPI_AFLOAT
 
#endif

#ifdef STDC_HEADERS
int HPLAI_bcast_blong
(
   HPLAI_T_panel                * PANEL,
   int                        * IFLAG
)
#else
int HPLAI_bcast_blong( PANEL, IFLAG )
   HPLAI_T_panel                * PANEL;
   int                        * IFLAG;
#endif
{ 
/*
 * .. Local Variables ..
 */
   MPI_Comm                   comm;
   int                        COUNT, count, dummy=0, ierr=MPI_SUCCESS,
                              ibuf, ibufR, ibufS, indx, ip2, k, l, lbuf,
                              lbufR, lbufS, mask, msgid, mydist, mydist2,
                              next, npm1, partner, prev, rank, root, size;
/* ..
 * .. Executable Statements ..
 */
   if( PANEL == NULL ) { *IFLAG = HPL_SUCCESS; return( HPL_SUCCESS ); }
   if( ( size = PANEL->grid->npcol ) <= 1 )
   {                     *IFLAG = HPL_SUCCESS; return( HPL_SUCCESS ); }
/*
 * Cast phase:  If I am the root process, start spreading the panel.  If
 * I am not the root process,  test  for  message receive completion. If
 * the message  is there,  then receive it,  and  keep  spreading  in  a
 * blocking fashion this time.  Otherwise,  inform  the caller  that the
 * panel has still not been received. 
 */
   comm    = PANEL->grid->row_comm;  rank  = PANEL->grid->mycol;
   mask    = PANEL->grid->col_mask;  ip2   = PANEL->grid->col_ip2m1;
   root    = PANEL->pcol;            msgid = PANEL->msgid;
   COUNT   = PANEL->len;             npm1  = size - 1;
   mydist2 = ( mydist = MModSub( rank, root, size ) ); indx = ip2;
   count   = COUNT / size; count = Mmax( count, 1 );
/*
 * Spread the panel across process columns
 */
   do
   {
      mask ^= ip2;
 
      if( ( mydist & mask ) == 0 )
      {
         lbuf = COUNT - ( ibuf = indx * count );
         if( indx + ip2 < size ) { l = ip2 * count; lbuf = Mmin( lbuf, l ); }
 
         partner = mydist ^ ip2;
 
         if( ( mydist & ip2 ) != 0 )
         {
            partner = MModAdd( root, partner, size );
/*
 * This probing mechanism causes problems when lookhead is on. Too many
 * messages are exchanged  in this virtual topology  causing  a hang on 
 * some machines. It is currently disabled until a better understanding
 * is acquired.
 */
#if 0
            ierr = MPI_Iprobe( partner, msgid, comm, &go, &PANEL->status[0] );
            if( ierr == MPI_SUCCESS )
            {        /* if panel is not here, return and keep testing */
               if( go == 0 )
               { *IFLAG = HPLAI_KEEP_TESTING; return( HPLAI_KEEP_TESTING ); }
            }
#endif
            if( lbuf > 0 )
            {
#ifdef HPL_USE_MPI_DATATYPE
               if( ierr == MPI_SUCCESS )
                  ierr =   HPLAI_packL( PANEL, ibuf, lbuf, I_RECV );
#endif
               if( ierr == MPI_SUCCESS )
                  ierr =   MPI_Recv( _M_BUFF_R, _M_COUNT_R, _M_TYPE_R,
                                     partner, msgid, comm, &PANEL->status[0] );
#ifdef HPL_USE_MPI_DATATYPE
               if( ierr == MPI_SUCCESS )
                  ierr =   MPI_Type_free( &PANEL->dtypes[I_RECV] );
#endif
            }
            else       /* Recv message of length zero to enable probe */
            {
               if( ierr == MPI_SUCCESS )
                  ierr =   MPI_Recv( (void *)(&dummy), 0, MPI_BYTE, partner,
                                     msgid, comm, &PANEL->status[0] );
            }
         }
         else if( partner < size )
         {
            partner = MModAdd( root, partner, size );
 
            if( lbuf > 0 )
            {
#ifdef HPL_USE_MPI_DATATYPE
               if( ierr == MPI_SUCCESS )
                  ierr =   HPLAI_packL( PANEL, ibuf, lbuf, I_SEND );
#endif
               if( ierr == MPI_SUCCESS )
                  ierr =   MPI_Ssend( _M_BUFF_S, _M_COUNT_S, _M_TYPE_S,
                                      partner, msgid, comm );
#ifdef HPL_USE_MPI_DATATYPE
               if( ierr == MPI_SUCCESS )
                  ierr =   MPI_Type_free( &PANEL->dtypes[I_SEND] );
#endif
            }
            else       /* Send message of length zero to enable probe */
            {
               if( ierr == MPI_SUCCESS )
                  ierr =   MPI_Ssend( (void *)(&dummy), 0, MPI_BYTE,
                                      partner, msgid, comm );
            }
         }
      }
 
      if( mydist2 < ip2 ) {  ip2 >>= 1; indx -= ip2; }
      else { mydist2 -= ip2; ip2 >>= 1; indx += ip2; }
 
   } while( ip2 > 0 );
/*
 * Roll the pieces
 */
   prev = MModSub1( rank, size ); next = MModAdd1( rank, size );

   for( k = 0; k < npm1; k++ )
   {
      l = ( k >> 1 ); 
/*
 * Who is sending to who and how much
 */
      if( ( ( mydist + k ) & 1 ) != 0 )
      {
         ibufS = ( indx = MModAdd( mydist, l,   size ) ) * count;
         lbufS = ( indx == npm1 ? COUNT : ibufS + count );
         lbufS = Mmin( COUNT, lbufS ) - ibufS; lbufS = Mmax( 0, lbufS );

         ibufR = ( indx = MModSub( mydist, l+1, size ) ) * count;
         lbufR = ( indx == npm1 ? COUNT : ibufR + count );
         lbufR = Mmin( COUNT, lbufR ) - ibufR; lbufR = Mmax( 0, lbufR );

         partner = prev;
      }
      else
      {
         ibufS = ( indx = MModSub( mydist, l,   size ) ) * count;
         lbufS = ( indx == npm1 ? COUNT : ibufS + count );
         lbufS = Mmin( COUNT, lbufS ) - ibufS; lbufS = Mmax( 0, lbufS );

         ibufR = ( indx = MModAdd( mydist, l+1, size ) ) * count;
         lbufR = ( indx == npm1 ? COUNT : ibufR + count );
         lbufR = Mmin( COUNT, lbufR ) - ibufR; lbufR = Mmax( 0, lbufR );

         partner = next;
      }
/*
 * Exchange the messages
 */
      if( lbufS > 0 )
      {
#ifdef HPL_USE_MPI_DATATYPE
         if( ierr == MPI_SUCCESS )
            ierr =   HPLAI_packL( PANEL, ibufS, lbufS, I_SEND );
#endif
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Issend( _M_ROLL_BUFF_S, _M_ROLL_COUNT_S,
                                 _M_ROLL_TYPE_S, partner, msgid, comm,
                                 &PANEL->request[0] );
      }
      else
      {
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Issend( (void *)(&dummy), 0, MPI_BYTE, partner,
                                 msgid, comm, &PANEL->request[0] );
      }

      if(  lbufR > 0 )
      {
#ifdef HPL_USE_MPI_DATATYPE
         if( ierr == MPI_SUCCESS )
            ierr =   HPLAI_packL( PANEL, ibufR, lbufR, I_RECV );
#endif
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Recv( _M_ROLL_BUFF_R, _M_ROLL_COUNT_R,
                               _M_ROLL_TYPE_R, partner, msgid, comm,
                               &PANEL->status[0] );
#ifdef HPL_USE_MPI_DATATYPE
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Type_free( &PANEL->dtypes[I_RECV] );
#endif
      }
      else
      {
         if( ierr == MPI_SUCCESS )
            ierr =   MPI_Recv( (void *)(&dummy), 0, MPI_BYTE, partner,
                               msgid, comm, &PANEL->status[0] );
      }

      if( ierr == MPI_SUCCESS )
         ierr =   MPI_Wait ( &PANEL->request[0], &PANEL->status[0] );
#ifdef HPL_USE_MPI_DATATYPE
      if( ( lbufS > 0 ) && ( ierr == MPI_SUCCESS ) )
         ierr =   MPI_Type_free( &PANEL->dtypes[I_SEND] );
#endif
   }
/*
 * If the message was received and being forwarded,  return HPL_SUCCESS.
 * If an error occured in an MPI call, return HPL_FAILURE.
 */
   *IFLAG = ( ierr == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE );

   return( *IFLAG );
}

#ifdef STDC_HEADERS
int HPLAI_bwait_blong
(
   HPLAI_T_panel *              PANEL
)
#else
int HPLAI_bwait_blong( PANEL )
   HPLAI_T_panel *              PANEL;
#endif
{
/* ..
 * .. Executable Statements ..
 */
   if( PANEL == NULL )           { return( HPL_SUCCESS ); }
   if( PANEL->grid->npcol <= 1 ) { return( HPL_SUCCESS ); }

   return( HPL_SUCCESS );
}

#ifdef __cplusplus
}
#endif
