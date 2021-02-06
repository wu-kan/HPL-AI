#ifndef HPLAI_PFACT_H
#define HPLAI_PFACT_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_pfact.h"
#include "hplai_misc.h"
#include "hplai_blas.h"

#include "hplai_panel.h"
/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */
typedef void (*HPLAI_T_PFA_FUN)
(  HPLAI_T_panel *,   const int,       const int,       const int,
   double * );
typedef void (*HPLAI_T_RFA_FUN)
(  HPLAI_T_panel *,   const int,       const int,       const int,
   double * );
typedef void (*HPLAI_T_UPD_FUN)
(  HPLAI_T_panel *,   int *,           HPLAI_T_panel *,   const int ); 

#endif
/*
 * End of hplai_pfact.h
 */
