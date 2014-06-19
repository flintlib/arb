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

#ifndef ACB_CALC_H
#define ACB_CALC_H

#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "arb_calc.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*acb_calc_func_t)(acb_ptr out,
    const acb_t inp, void * param, long order, long prec);

/* Bounds */

void acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func,
    void * param, const acb_t x, const arb_t radius,
    long maxdepth, long prec);

/* Integration */

int acb_calc_integrate_taylor(acb_t res,
    acb_calc_func_t func, void * param,
    const acb_t a, const acb_t b,
    const arf_t inner_radius,
    const arf_t outer_radius,
    long accuracy_goal, long prec);

#ifdef __cplusplus
}
#endif

#endif

