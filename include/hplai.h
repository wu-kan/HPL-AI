#ifndef HPLAI_H
#define HPLAI_H

#ifdef __cplusplus
extern "C"
{
#endif
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl.h"

/*
 * ---------------------------------------------------------------------
 * typedef definitions
 * ---------------------------------------------------------------------
 */
#ifndef HPLAI_T_AFLOAT
#define HPLAI_T_AFLOAT float
#endif

#include "hplai_blas.h"
#include "hplai_comm.h"
#include "hplai_grid.h"
#include "hplai_panel.h"
#include "hplai_pgesv.h"

#include "hplai_pmatgen.h"
#include "hplai_ptest.h"

#ifdef __cplusplus
}
#endif

#endif
/*
 * End of hplai.h
 */
