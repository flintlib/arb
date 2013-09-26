/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#ifndef FMPRB_CALC_H
#define FMPRB_CALC_H

#include "fmprb.h"
#include "fmprb_poly.h"
#include "fmprb_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

extern TLS_PREFIX int fmprb_calc_verbose;

typedef int (*fmprb_calc_func_t)(fmprb_ptr out,
    const fmprb_t inp, void * param, long order, long prec);

#define FMPRB_CALC_SUCCESS 0
#define FMPRB_CALC_IMPRECISE_INPUT 1
#define FMPRB_CALC_NO_CONVERGENCE 2

/* Root-finding */

long fmprb_calc_isolate_roots(fmprb_ptr * blocks, int ** flags,
    fmprb_calc_func_t func, void * param,
    const fmprb_t block, long maxdepth, long maxeval, long maxfound,
    long prec);

void fmprb_calc_newton_conv_factor(fmpr_t conv_factor,
    fmprb_calc_func_t func, void * param, const fmprb_t conv_region, long prec);

int fmprb_calc_newton_step(fmprb_t xnew, fmprb_calc_func_t func, void * param,
    const fmprb_t x, const fmprb_t conv_region, const fmpr_t conv_factor, long prec);

int fmprb_calc_refine_root_newton(fmprb_t r, fmprb_calc_func_t func,
    void * param, const fmprb_t start, const fmprb_t conv_region,
    const fmpr_t conv_factor, long eval_extra_prec, long prec);

#ifdef __cplusplus
}
#endif

#endif

