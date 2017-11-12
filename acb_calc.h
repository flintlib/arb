/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ACB_CALC_H
#define ACB_CALC_H

#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "arb_calc.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ACB_CALC_VERBOSE  1
#define ACB_CALC_VERY_VERBOSE 3

#define ACB_CALC_INTEGRATE_HEAP 8

typedef int (*acb_calc_func_t)(acb_ptr out,
    const acb_t inp, void * param, slong order, slong prec);

/* Bounds */

void acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func,
    void * param, const acb_t x, const arb_t radius,
    slong maxdepth, slong prec);

/* Integration */

int acb_calc_integrate_taylor(acb_t res,
    acb_calc_func_t func, void * param,
    const acb_t a, const acb_t b,
    const arf_t inner_radius,
    const arf_t outer_radius,
    slong accuracy_goal, slong prec);

void
acb_calc_integrate(acb_t res, acb_calc_func_t f, void * param,
    const acb_t a, const acb_t b,
    slong goal, const mag_t tol,
    slong deg_limit, slong eval_limit, slong depth_limit,
    int flags,
    slong prec);

slong
acb_calc_integrate_gl_auto_deg(acb_t res, acb_calc_func_t f, void * param,
        const acb_t a, const acb_t b, const mag_t tol,
        slong deg_limit, int flags, slong prec);

#ifdef __cplusplus
}
#endif

#endif

