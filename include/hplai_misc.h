#ifndef HPLAI_MISC_H
#define HPLAI_MISC_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */

#include <math.h>
#include <float.h>

/*
 * ---------------------------------------------------------------------
 * #define macros definitions
 * ---------------------------------------------------------------------
 */
#ifndef HPLAI_T_AFLOAT
#define HPLAI_T_AFLOAT float
#endif

#define HPLAI_rone ((HPLAI_T_AFLOAT)HPL_rone)
#define HPLAI_rtwo ((HPLAI_T_AFLOAT)HPL_rtwo)
#define HPLAI_rzero ((HPLAI_T_AFLOAT)HPL_rzero)

#define HPLAI_PTR HPL_PTR

extern MPI_Datatype HPLAI_MPI_AFLOAT;

void HPLAI_init
    STDC_ARGS((
        const int,
        const int));

void HPLAI_finalize();

#endif
/*
 * End of hpl_misc.h
 */
