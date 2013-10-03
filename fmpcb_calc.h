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

#ifndef FMPCB_CALC_H
#define FMPCB_CALC_H

#include "fmpcb.h"
#include "fmpcb_poly.h"
#include "fmpcb_mat.h"
#include "fmprb_calc.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*fmpcb_calc_func_t)(fmpcb_ptr out,
    const fmpcb_t inp, void * param, long order, long prec);

/* Bounds */

void fmpcb_calc_cauchy_bound(fmprb_t bound, fmpcb_calc_func_t func,
    void * param, const fmpcb_t x, const fmprb_t radius,
    long maxdepth, long prec);

/* Integration */

int fmpcb_calc_integrate_taylor(fmpcb_t res,
    fmpcb_calc_func_t func, void * param,
    const fmpcb_t a, const fmpcb_t b,
    const fmpr_t inner_radius,
    const fmpr_t outer_radius,
    long accuracy_goal, long prec);

#ifdef __cplusplus
}
#endif

#endif

