#ifndef HPLAI_BLAS_H
#define HPLAI_BLAS_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_blas.h"
#include "hplai_misc.h"

#define HPLAI_ORDER HPL_ORDER
#define HPLAI_RowMajor HplRowMajor
#define HPLAI_ColumnMajor HplColumnMajor

#define HPLAI_TRANS HPL_TRANS
#define HPLAI_NoTrans HplNoTrans
#define HPLAI_Trans HplTrans
#define HPLAI_ConjTrans HplConjTrans

#define HPLAI_UPLO HPL_UPLO
#define HPLAI_Upper HplUpper
#define HPLAI_Lower HplLower

#define HPLAI_DIAG HPL_DIAG
#define HPLAI_NonUnit HplNonUnit
#define HPLAI_Unit HplUnit

#define HPLAI_SIDE HPL_SIDE
#define HPLAI_Left HplLeft
#define HPLAI_Right HplRight

int HPLAI_iaamax
    STDC_ARGS((
        const int,
        const HPLAI_T_AFLOAT *,
        const int));
void HPLAI_aaxpy
    STDC_ARGS((
        const int,
        const HPLAI_T_AFLOAT,
        const HPLAI_T_AFLOAT *,
        const int,
        HPLAI_T_AFLOAT *,
        const int));
void HPLAI_acopy
    STDC_ARGS((
        const int,
        const HPLAI_T_AFLOAT *,
        const int,
        HPLAI_T_AFLOAT *,
        const int));
void HPLAI_ascal
    STDC_ARGS((
        const int,
        const HPLAI_T_AFLOAT,
        HPLAI_T_AFLOAT *,
        const int));
void HPLAI_agemv
    STDC_ARGS((
        const enum HPL_ORDER,
        const enum HPL_TRANS,
        const int,
        const int,
        const HPLAI_T_AFLOAT,
        const HPLAI_T_AFLOAT *,
        const int,
        const HPLAI_T_AFLOAT *,
        const int,
        const HPLAI_T_AFLOAT,
        HPLAI_T_AFLOAT *,
        const int));
void HPLAI_ager
    STDC_ARGS((
        const enum HPL_ORDER,
        const int,
        const int,
        const HPLAI_T_AFLOAT,
        const HPLAI_T_AFLOAT *,
        const int,
        HPLAI_T_AFLOAT *,
        const int,
        HPLAI_T_AFLOAT *,
        const int));
void HPLAI_atrsv
    STDC_ARGS((
        const enum HPL_ORDER,
        const enum HPL_UPLO,
        const enum HPL_TRANS,
        const enum HPL_DIAG,
        const int,
        const HPLAI_T_AFLOAT *,
        const int,
        HPLAI_T_AFLOAT *,
        const int));
void HPLAI_agemm
    STDC_ARGS((
        const enum HPL_ORDER,
        const enum HPL_TRANS,
        const enum HPL_TRANS,
        const int,
        const int,
        const int,
        const HPLAI_T_AFLOAT,
        const HPLAI_T_AFLOAT *,
        const int,
        const HPLAI_T_AFLOAT *,
        const int,
        const HPLAI_T_AFLOAT,
        HPLAI_T_AFLOAT *,
        const int));
void HPLAI_atrsm
    STDC_ARGS((
        const enum HPL_ORDER,
        const enum HPL_SIDE,
        const enum HPL_UPLO,
        const enum HPL_TRANS,
        const enum HPL_DIAG,
        const int,
        const int,
        const HPLAI_T_AFLOAT,
        const HPLAI_T_AFLOAT *,
        const int,
        HPLAI_T_AFLOAT *,
        const int));

#endif
/*
 * hplai_blas.h
 */
