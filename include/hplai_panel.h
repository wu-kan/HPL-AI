#ifndef HPLAI_PANEL_H
#define HPLAI_PANEL_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_panel.h"
#include "hplai_misc.h"
#include "hplai_grid.h"
#include "hplai_pfact.h"
/*
 * ---------------------------------------------------------------------
 * panel function prototypes
 * ---------------------------------------------------------------------
 */

void                             HPLAI_papanel_new
STDC_ARGS( (
   HPL_T_grid *,
   HPLAI_T_palg *,
   const int,
   const int,
   const int,
   HPLAI_T_pmat *,
   const int,
   const int,
   const int,
   HPLAI_T_panel * *
) );
void                             HPLAI_papanel_init
STDC_ARGS( (
   HPL_T_grid *,
   HPLAI_T_palg *,
   const int,
   const int,
   const int,
   HPLAI_T_pmat *,
   const int,
   const int,
   const int,
   HPLAI_T_panel *
) );
int                              HPLAI_papanel_disp
STDC_ARGS( (
   HPLAI_T_panel * *
) );
int                              HPLAI_papanel_free
STDC_ARGS( (
   HPLAI_T_panel *
) );

#endif
/*
 * End of hplai_panel.h
 */
