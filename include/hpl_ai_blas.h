#ifndef HPL_AI_BLAS_H
#define HPL_AI_BLAS_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_ai.h"

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

void HPL_AI_blas_init
    STDC_ARGS((
        const int,
        const int));

void HPL_AI_blas_finalize();

int HPL_iaamax
    STDC_ARGS((
        const int,
        const HPL_AI_T_afloat *,
        const int));
void HPL_AI_aaxpy
    STDC_ARGS((
        const int,
        const HPL_AI_T_afloat,
        const HPL_AI_T_afloat *,
        const int,
        HPL_AI_T_afloat *,
        const int));
void HPL_AI_acopy
    STDC_ARGS((
        const int,
        const HPL_AI_T_afloat *,
        const int,
        HPL_AI_T_afloat *,
        const int));
void HPL_AI_ascal
    STDC_ARGS((
        const int,
        const HPL_AI_T_afloat,
        HPL_AI_T_afloat *,
        const int));
void HPL_AI_aswap
    STDC_ARGS((
        const int,
        HPL_AI_T_afloat *,
        const int,
        HPL_AI_T_afloat *,
        const int));
void HPL_AI_agemv
    STDC_ARGS((
        const enum HPL_ORDER,
        const enum HPL_TRANS,
        const int,
        const int,
        const HPL_AI_T_afloat,
        const HPL_AI_T_afloat *,
        const int,
        const HPL_AI_T_afloat *,
        const int,
        const HPL_AI_T_afloat,
        HPL_AI_T_afloat *,
        const int));
void HPL_AI_ager
    STDC_ARGS((
        const enum HPL_ORDER,
        const int,
        const int,
        const HPL_AI_T_afloat,
        const HPL_AI_T_afloat *,
        const int,
        HPL_AI_T_afloat *,
        const int,
        HPL_AI_T_afloat *,
        const int));
void HPL_AI_atrsv
    STDC_ARGS((
        const enum HPL_ORDER,
        const enum HPL_UPLO,
        const enum HPL_TRANS,
        const enum HPL_DIAG,
        const int,
        const HPL_AI_T_afloat *,
        const int,
        HPL_AI_T_afloat *,
        const int));
void HPL_AI_agemm
    STDC_ARGS((
        const enum HPL_ORDER,
        const enum HPL_TRANS,
        const enum HPL_TRANS,
        const int,
        const int,
        const int,
        const HPL_AI_T_afloat,
        const HPL_AI_T_afloat *,
        const int,
        const HPL_AI_T_afloat *,
        const int,
        const HPL_AI_T_afloat,
        HPL_AI_T_afloat *,
        const int));
void HPL_AI_atrsm
    STDC_ARGS((
        const enum HPL_ORDER,
        const enum HPL_SIDE,
        const enum HPL_UPLO,
        const enum HPL_TRANS,
        const enum HPL_DIAG,
        const int,
        const int,
        const HPL_AI_T_afloat,
        const HPL_AI_T_afloat *,
        const int,
        HPL_AI_T_afloat *,
        const int));

#endif
/*
 * hpl_ai_blas.h
 */
