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

#include "fmprb_calc.h"

int fmprb_calc_newton_step(fmprb_t xnew, fmprb_calc_func_t func,
    void * param, const fmprb_t x, const fmprb_t conv_region,
    const fmpr_t conv_factor, long prec)
{
    fmpr_t err;
    fmprb_t t;
    fmprb_struct u[2];
    int result;

    fmpr_init(err);
    fmprb_init(t);
    fmprb_init(u + 0);
    fmprb_init(u + 1);

    fmpr_mul(err, fmprb_radref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_mul(err, err, conv_factor, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmpr_set(fmprb_midref(t), fmprb_midref(x));
    fmpr_zero(fmprb_radref(t));

    func(u, t, param, 2, prec);

    fmprb_div(u, u, u + 1, prec);
    fmprb_sub(u, t, u, prec);

    fmprb_add_error_fmpr(u, err);

    if (fmprb_contains(conv_region, u) &&
        (fmpr_cmp(fmprb_radref(u), fmprb_radref(x)) < 0))
    {
        fmprb_swap(xnew, u);
        result = FMPRB_CALC_SUCCESS;
    }
    else
    {
        fmprb_set(xnew, x);
        result = FMPRB_CALC_NO_CONVERGENCE;
    }

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(u + 1);
    fmpr_clear(err);

    return result;
}

