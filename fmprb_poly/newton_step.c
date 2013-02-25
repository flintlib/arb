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

#include "fmprb_poly.h"

int
_fmprb_poly_newton_step(fmprb_t xnew, const fmprb_struct * poly, long len,
    const fmprb_t x,
    const fmprb_t convergence_interval,
    const fmpr_t convergence_factor, long prec)
{
    fmpr_t err;
    fmprb_t t, u, v;
    int result;

    fmpr_init(err);
    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    fmpr_mul(err, fmprb_radref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_mul(err, err, convergence_factor, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmpr_set(fmprb_midref(t), fmprb_midref(x));
    fmpr_zero(fmprb_radref(t));

    _fmprb_poly_evaluate2(u, v, poly, len, t, prec);

    fmprb_div(u, u, v, prec);
    fmprb_sub(u, t, u, prec);

    fmprb_add_error_fmpr(u, err);

    if (fmprb_contains(convergence_interval, u) &&
        (fmpr_cmp(fmprb_radref(u), fmprb_radref(x)) < 0))
    {
        fmprb_swap(xnew, u);
        result = 1;
    }
    else
    {
        fmprb_set(xnew, x);
        result = 0;
    }

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
    fmpr_clear(err);

    return result;
}

