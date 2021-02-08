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
#define HPLAI_TRANS HPL_TRANS
#define HPLAI_UPLO HPL_UPLO
#define HPLAI_DIAG HPL_DIAG
#define HPLAI_SIDE HPL_SIDE

#if (HPL_CALL_CBLAS)

//float---
CBLAS_INDEX       cblas_isamax
STDC_ARGS(
(  const int,       const float *,  const int ) );
void              cblas_sswap
STDC_ARGS(
(  const int,       float *,        const int,       float *,
   const int ) );
void              cblas_scopy
STDC_ARGS(
(  const int,       const float *,  const int,       float *,
   const int ) );
void              cblas_saxpy
STDC_ARGS(
(  const int,       const float,    const float *,  const int,
   float *,        const int ) );
void              cblas_sscal
STDC_ARGS(
(  const int,       const float,    float *,        const int ) );

void              cblas_sgemv
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_TRANSPOSE,
   const int,       const int,       const float,    const float *,
   const int,       const float *,  const int,       const float,
   float *,        const int ) );

void              cblas_sger
STDC_ARGS(
(  const enum CBLAS_ORDER,           const int,       const int,
   const float,    const float *,  const int,       const float *,
   const int,       float *,        const int ) );
void              cblas_strsv
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_UPLO,
   const enum CBLAS_TRANSPOSE,       const enum CBLAS_DIAG,
   const int,       const float *,  const int,       float *,
   const int ) );

void              cblas_sgemm
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_TRANSPOSE,
   const enum CBLAS_TRANSPOSE,       const int,       const int,
   const int,       const float,    const float *,  const int,
   const float *,  const int,       const float,    float *,
   const int ) );
void              cblas_strsm
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_SIDE,
   const enum CBLAS_UPLO,            const enum CBLAS_TRANSPOSE,
   const enum CBLAS_DIAG,            const int,       const int,
   const float,    const float *,  const int,       float *,
   const int ) );
#endif

void HPLAI_blas_init
    STDC_ARGS((
        const int,
        const int));

void HPLAI_blas_finalize();

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
