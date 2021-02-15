/*
 * Include files
 */
#include "hplai.h"
/*
 * Define default value for unrolling factors
 * #ifndef HPL_LACPY_M_DEPTH
 * #define    HPL_LACPY_M_DEPTH       32
 * #define    HPL_LACPY_LOG2_M_DEPTH   5
 * #endif
 * #ifndef HPL_LACPY_N_DEPTH
 * #define    HPL_LACPY_N_DEPTH        4
 * #define    HPL_LACPY_LOG2_N_DEPTH   2
 * #endif
 */
#ifndef HPL_LACPY_M_DEPTH
#define    HPL_LACPY_M_DEPTH        4
#define    HPL_LACPY_LOG2_M_DEPTH   2
#endif
#ifndef HPL_LACPY_N_DEPTH
#define    HPL_LACPY_N_DEPTH        2
#define    HPL_LACPY_LOG2_N_DEPTH   1
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
   void HPLAI_alacpy(
       const int M,
       const int N,
       const HPLAI_T_AFLOAT *A,
       const int LDA,
       HPLAI_T_AFLOAT *B,
       const int LDB)
#else
void HPLAI_alacpy(M, N, A, LDA, B, LDB)
    const int M;
const int N;
const HPLAI_T_AFLOAT *A;
const int LDA;
HPLAI_T_AFLOAT *B;
const int LDB;
#endif
   {
/*
 * .. Local Variables ..
 */
#ifdef HPL_LACPY_USE_COPY
   register int               j;
#else
#if   ( HPL_LACPY_N_DEPTH ==  1 )
   const HPLAI_T_AFLOAT               * A0 = A;
   HPLAI_T_AFLOAT                     * B0 = B;
#elif ( HPL_LACPY_N_DEPTH ==  2 )
   const HPLAI_T_AFLOAT               * A0 = A,              * A1 = A +     LDA;
   HPLAI_T_AFLOAT                     * B0 = B,              * B1 = B +     LDB;
#elif ( HPL_LACPY_N_DEPTH ==  4 )
   const HPLAI_T_AFLOAT               * A0 = A,              * A1 = A +     LDA,
                              * A2 = A + (LDA << 1), * A3 = A + 3 * LDA;
   HPLAI_T_AFLOAT                     * B0 = B,              * B1 = B +     LDB,
                              * B2 = B + (LDB << 1), * B3 = B + 3 * LDB;
#endif
   const int                  incA = ( (unsigned int)(LDA) <<
                                       HPL_LACPY_LOG2_N_DEPTH ) - M,
                              incB = ( (unsigned int)(LDB) <<
                                       HPL_LACPY_LOG2_N_DEPTH ) - M,
                              incA0 = (unsigned int)(LDA) - M,
                              incB0 = (unsigned int)(LDB) - M;
   int                        mu, nu;
   register int               i, j;
#endif
/* ..
 * .. Executable Statements ..
 */
   if( ( M <= 0 ) || ( N <= 0 ) ) return;

#ifdef HPL_LACPY_USE_COPY
   for( j = 0; j < N; j++, A0 += LDA, B0 += LDB ) HPLAI_acopy( M, A0, 1, B0, 1 );
#else
   mu = (int)( ( (unsigned int)(M) >> HPL_LACPY_LOG2_M_DEPTH ) <<
                                      HPL_LACPY_LOG2_M_DEPTH );
   nu = (int)( ( (unsigned int)(N) >> HPL_LACPY_LOG2_N_DEPTH ) <<
                                      HPL_LACPY_LOG2_N_DEPTH );

   for( j = 0; j < nu; j += HPL_LACPY_N_DEPTH )
   {
      for( i = 0; i < mu; i += HPL_LACPY_M_DEPTH )
      {
#if   ( HPL_LACPY_N_DEPTH ==  1 )
         B0[ 0] = A0[ 0];
#elif ( HPL_LACPY_N_DEPTH ==  2 )
         B0[ 0] = A0[ 0]; B1[ 0] = A1[ 0];
#elif ( HPL_LACPY_N_DEPTH ==  4 )
         B0[ 0] = A0[ 0]; B1[ 0] = A1[ 0]; B2[ 0] = A2[ 0]; B3[ 0] = A3[ 0];
#endif

#if ( HPL_LACPY_M_DEPTH >  1 )

#if   ( HPL_LACPY_N_DEPTH ==  1 )
         B0[ 1] = A0[ 1];
#elif ( HPL_LACPY_N_DEPTH ==  2 )
         B0[ 1] = A0[ 1]; B1[ 1] = A1[ 1];
#elif ( HPL_LACPY_N_DEPTH ==  4 )
         B0[ 1] = A0[ 1]; B1[ 1] = A1[ 1]; B2[ 1] = A2[ 1]; B3[ 1] = A3[ 1];
#endif

#endif
#if ( HPL_LACPY_M_DEPTH >  2 )

#if   ( HPL_LACPY_N_DEPTH ==  1 )
         B0[ 2] = A0[ 2]; B0[ 3] = A0[ 3];
#elif ( HPL_LACPY_N_DEPTH ==  2 )
         B0[ 2] = A0[ 2]; B1[ 2] = A1[ 2]; B0[ 3] = A0[ 3]; B1[ 3] = A1[ 3];
#elif ( HPL_LACPY_N_DEPTH ==  4 )
         B0[ 2] = A0[ 2]; B1[ 2] = A1[ 2]; B2[ 2] = A2[ 2]; B3[ 2] = A3[ 2];
         B0[ 3] = A0[ 3]; B1[ 3] = A1[ 3]; B2[ 3] = A2[ 3]; B3[ 3] = A3[ 3];
#endif

#endif
#if ( HPL_LACPY_M_DEPTH >  4 )

#if   ( HPL_LACPY_N_DEPTH ==  1 )
         B0[ 4] = A0[ 4]; B0[ 5] = A0[ 5]; B0[ 6] = A0[ 6]; B0[ 7] = A0[ 7];
#elif ( HPL_LACPY_N_DEPTH ==  2 )
         B0[ 4] = A0[ 4]; B1[ 4] = A1[ 4]; B0[ 5] = A0[ 5]; B1[ 5] = A1[ 5];
         B0[ 6] = A0[ 6]; B1[ 6] = A1[ 6]; B0[ 7] = A0[ 7]; B1[ 7] = A1[ 7];
#elif ( HPL_LACPY_N_DEPTH ==  4 )
         B0[ 4] = A0[ 4]; B1[ 4] = A1[ 4]; B2[ 4] = A2[ 4]; B3[ 4] = A3[ 4];
         B0[ 5] = A0[ 5]; B1[ 5] = A1[ 5]; B2[ 5] = A2[ 5]; B3[ 5] = A3[ 5];
         B0[ 6] = A0[ 6]; B1[ 6] = A1[ 6]; B2[ 6] = A2[ 6]; B3[ 6] = A3[ 6];
         B0[ 7] = A0[ 7]; B1[ 7] = A1[ 7]; B2[ 7] = A2[ 7]; B3[ 7] = A3[ 7];
#endif

#endif
#if ( HPL_LACPY_M_DEPTH >  8 )

#if   ( HPL_LACPY_N_DEPTH ==  1 )
         B0[ 8] = A0[ 8]; B0[ 9] = A0[ 9]; B0[10] = A0[10]; B0[11] = A0[11];
         B0[12] = A0[12]; B0[13] = A0[13]; B0[14] = A0[14]; B0[15] = A0[15];
#elif ( HPL_LACPY_N_DEPTH ==  2 )
         B0[ 8] = A0[ 8]; B1[ 8] = A1[ 8]; B0[ 9] = A0[ 9]; B1[ 9] = A1[ 9];
         B0[10] = A0[10]; B1[10] = A1[10]; B0[11] = A0[11]; B1[11] = A1[11];
         B0[12] = A0[12]; B1[12] = A1[12]; B0[13] = A0[13]; B1[13] = A1[13];
         B0[14] = A0[14]; B1[14] = A1[14]; B0[15] = A0[15]; B1[15] = A1[15];
#elif ( HPL_LACPY_N_DEPTH ==  4 )
         B0[ 8] = A0[ 8]; B1[ 8] = A1[ 8]; B2[ 8] = A2[ 8]; B3[ 8] = A3[ 8];
         B0[ 9] = A0[ 9]; B1[ 9] = A1[ 9]; B2[ 9] = A2[ 9]; B3[ 9] = A3[ 9];
         B0[10] = A0[10]; B1[10] = A1[10]; B2[10] = A2[10]; B3[10] = A3[10];
         B0[11] = A0[11]; B1[11] = A1[11]; B2[11] = A2[11]; B3[11] = A3[11];
         B0[12] = A0[12]; B1[12] = A1[12]; B2[12] = A2[12]; B3[12] = A3[12];
         B0[13] = A0[13]; B1[13] = A1[13]; B2[13] = A2[13]; B3[13] = A3[13];
         B0[14] = A0[14]; B1[14] = A1[14]; B2[14] = A2[14]; B3[14] = A3[14];
         B0[15] = A0[15]; B1[15] = A1[15]; B2[15] = A2[15]; B3[15] = A3[15];
#endif

#endif
#if ( HPL_LACPY_M_DEPTH > 16 )

#if   ( HPL_LACPY_N_DEPTH ==  1 )
         B0[16] = A0[16]; B0[17] = A0[17]; B0[18] = A0[18]; B0[19] = A0[19];
         B0[20] = A0[20]; B0[21] = A0[21]; B0[22] = A0[22]; B0[23] = A0[23];
         B0[24] = A0[24]; B0[25] = A0[25]; B0[26] = A0[26]; B0[27] = A0[27];
         B0[28] = A0[28]; B0[29] = A0[29]; B0[30] = A0[30]; B0[31] = A0[31];
#elif ( HPL_LACPY_N_DEPTH ==  2 )
         B0[16] = A0[16]; B1[16] = A1[16]; B0[17] = A0[17]; B1[17] = A1[17];
         B0[18] = A0[18]; B1[18] = A1[18]; B0[19] = A0[19]; B1[19] = A1[19];
         B0[20] = A0[20]; B1[20] = A1[20]; B0[21] = A0[21]; B1[21] = A1[21];
         B0[22] = A0[22]; B1[22] = A1[22]; B0[23] = A0[23]; B1[23] = A1[23];
         B0[24] = A0[24]; B1[24] = A1[24]; B0[25] = A0[25]; B1[25] = A1[25];
         B0[26] = A0[26]; B1[26] = A1[26]; B0[27] = A0[27]; B1[27] = A1[27];
         B0[28] = A0[28]; B1[28] = A1[28]; B0[29] = A0[29]; B1[29] = A1[29];
         B0[30] = A0[30]; B1[30] = A1[30]; B0[31] = A0[31]; B1[31] = A1[31];
#elif ( HPL_LACPY_N_DEPTH ==  4 )
         B0[16] = A0[16]; B1[16] = A1[16]; B2[16] = A2[16]; B3[16] = A3[16];
         B0[17] = A0[17]; B1[17] = A1[17]; B2[17] = A2[17]; B3[17] = A3[17];
         B0[18] = A0[18]; B1[18] = A1[18]; B2[18] = A2[18]; B3[18] = A3[18];
         B0[19] = A0[19]; B1[19] = A1[19]; B2[19] = A2[19]; B3[19] = A3[19];
         B0[20] = A0[20]; B1[20] = A1[20]; B2[20] = A2[20]; B3[20] = A3[20];
         B0[21] = A0[21]; B1[21] = A1[21]; B2[21] = A2[21]; B3[21] = A3[21];
         B0[22] = A0[22]; B1[22] = A1[22]; B2[22] = A2[22]; B3[22] = A3[22];
         B0[23] = A0[23]; B1[23] = A1[23]; B2[23] = A2[23]; B3[23] = A3[23];
         B0[24] = A0[24]; B1[24] = A1[24]; B2[24] = A2[24]; B3[24] = A3[24];
         B0[25] = A0[25]; B1[25] = A1[25]; B2[25] = A2[25]; B3[25] = A3[25];
         B0[26] = A0[26]; B1[26] = A1[26]; B2[26] = A2[26]; B3[26] = A3[26];
         B0[27] = A0[27]; B1[27] = A1[27]; B2[27] = A2[27]; B3[27] = A3[27];
         B0[28] = A0[28]; B1[28] = A1[28]; B2[28] = A2[28]; B3[28] = A3[28];
         B0[29] = A0[29]; B1[29] = A1[29]; B2[29] = A2[29]; B3[29] = A3[29];
         B0[30] = A0[30]; B1[30] = A1[30]; B2[30] = A2[30]; B3[30] = A3[30];
         B0[31] = A0[31]; B1[31] = A1[31]; B2[31] = A2[31]; B3[31] = A3[31];
#endif

#endif

#if   ( HPL_LACPY_N_DEPTH ==  1 )
         A0 += HPL_LACPY_M_DEPTH; B0 += HPL_LACPY_M_DEPTH;
#elif ( HPL_LACPY_N_DEPTH ==  2 )
         A0 += HPL_LACPY_M_DEPTH; B0 += HPL_LACPY_M_DEPTH;
         A1 += HPL_LACPY_M_DEPTH; B1 += HPL_LACPY_M_DEPTH;
#elif ( HPL_LACPY_N_DEPTH ==  4 )
         A0 += HPL_LACPY_M_DEPTH; B0 += HPL_LACPY_M_DEPTH;
         A1 += HPL_LACPY_M_DEPTH; B1 += HPL_LACPY_M_DEPTH;
         A2 += HPL_LACPY_M_DEPTH; B2 += HPL_LACPY_M_DEPTH;
         A3 += HPL_LACPY_M_DEPTH; B3 += HPL_LACPY_M_DEPTH;
#endif
      }

      for( i = mu; i < M; i++ )
      {
#if   ( HPL_LACPY_N_DEPTH ==  1 )
         *B0 = *A0; B0++; A0++;
#elif ( HPL_LACPY_N_DEPTH ==  2 )
         *B0 = *A0; B0++; A0++; *B1 = *A1; B1++; A1++;
#elif ( HPL_LACPY_N_DEPTH ==  4 )
         *B0 = *A0; B0++; A0++; *B1 = *A1; B1++; A1++;
         *B2 = *A2; B2++; A2++; *B3 = *A3; B3++; A3++;
#endif
      }

#if   ( HPL_LACPY_N_DEPTH ==  1 )
      A0 += incA; B0 += incB;
#elif ( HPL_LACPY_N_DEPTH ==  2 )
      A0 += incA; B0 += incB; A1 += incA; B1 += incB;
#elif ( HPL_LACPY_N_DEPTH ==  4 )
      A0 += incA; B0 += incB; A1 += incA; B1 += incB;
      A2 += incA; B2 += incB; A3 += incA; B3 += incB;
#endif
   }

   for( j = nu; j < N; j++, B0 += incB0, A0 += incA0 )
   {
      for( i = 0; i < mu; i += HPL_LACPY_M_DEPTH,
           B0 += HPL_LACPY_M_DEPTH, A0 += HPL_LACPY_M_DEPTH )
      {
         B0[ 0] = A0[ 0];
#if ( HPL_LACPY_M_DEPTH >  1 )
         B0[ 1] = A0[ 1];
#endif
#if ( HPL_LACPY_M_DEPTH >  2 )
         B0[ 2] = A0[ 2]; B0[ 3] = A0[ 3];
#endif
#if ( HPL_LACPY_M_DEPTH >  4 )
         B0[ 4] = A0[ 4]; B0[ 5] = A0[ 5]; B0[ 6] = A0[ 6]; B0[ 7] = A0[ 7];
#endif
#if ( HPL_LACPY_M_DEPTH >  8 )
         B0[ 8] = A0[ 8]; B0[ 9] = A0[ 9]; B0[10] = A0[10]; B0[11] = A0[11];
         B0[12] = A0[12]; B0[13] = A0[13]; B0[14] = A0[14]; B0[15] = A0[15];
#endif
#if ( HPL_LACPY_M_DEPTH > 16 )
         B0[16] = A0[16]; B0[17] = A0[17]; B0[18] = A0[18]; B0[19] = A0[19];
         B0[20] = A0[20]; B0[21] = A0[21]; B0[22] = A0[22]; B0[23] = A0[23];
         B0[24] = A0[24]; B0[25] = A0[25]; B0[26] = A0[26]; B0[27] = A0[27];
         B0[28] = A0[28]; B0[29] = A0[29]; B0[30] = A0[30]; B0[31] = A0[31];
#endif
      }
      for( i = mu; i < M; i++, B0++, A0++ ) { *B0 = *A0; }
   }
#endif
   }

#ifdef __cplusplus
}
#endif
