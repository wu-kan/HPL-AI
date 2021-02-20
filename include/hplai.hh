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

#include "hplai_misc.h"

#include "hplai_auxil.h"
#include "hplai_gesv.h"

#include "hplai_pauxil.h"
#include "hplai_panel.h"
#include "hplai_pfact.h"

#include "hplai_pmatgen.h"
#include "hplai_ptest.h"

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

//use blaspp https://bitbucket.org/icl/blaspp/src/master/
#include "hplai_blas.hh"
#include "hplai_pgesv.hh"

#endif

#endif
/*
 * End of hplai.h
 */
