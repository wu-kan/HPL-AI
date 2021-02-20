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
   if( COUNT > 0 )
      return MPI_Bcast(BUFFER, COUNT, HPLAI_MPI_AFLOAT, ROOT, COMM);
/*
 * .. Local Variables ..
 */
   int                        hplerr=MPI_SUCCESS, ip2=1, kk, mask=1, 
                              mpierr, mydist, partner, rank, size, 
                              tag = MSGID_BEGIN_COLL;
   MPI_Status                 status;
/* ..
 * .. Executable Statements ..
 */
   if( COUNT <= 0 ) return( MPI_SUCCESS );
   mpierr = MPI_Comm_size( COMM, &size ); if( size <= 1 ) return( mpierr );
   mpierr = MPI_Comm_rank( COMM, &rank );

   kk = size - 1;
   while( kk > 1 ) { kk >>= 1; ip2 <<= 1; mask <<= 1; mask++; }
   mydist = MModSub( rank, ROOT, size );

   do
   {
      mask ^= ip2;
      if( ( mydist & mask ) == 0 )
      {
         partner = mydist ^ ip2;

         if( mydist & ip2 )
         {
            partner = MModAdd( ROOT, partner, size );
            mpierr  = MPI_Recv(  BUFFER, COUNT, HPLAI_MPI_AFLOAT,
                                 partner, tag, COMM, &status );
         }
         else if( partner < size )
         {
            partner = MModAdd( ROOT, partner, size );
            mpierr  = MPI_Send( BUFFER, COUNT, HPLAI_MPI_AFLOAT,
                                partner, tag, COMM );
         }
         if( mpierr != MPI_SUCCESS ) hplerr = mpierr;
      }
      ip2 >>= 1;
   } while( ip2 );

   return( hplerr );
}

#ifdef __cplusplus
}
#endif
