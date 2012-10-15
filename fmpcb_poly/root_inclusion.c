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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpcb_poly.h"

/*
Given any complex number m, and a polynomial f and its derivative f',
sets r to an interval centered on m that is guaranteed to contain
at least one root of f. Assumes len > 1.
*/
void
_fmpcb_poly_root_inclusion(fmpcb_t r, const fmpcb_t m,
    const fmpcb_struct * poly,
    const fmpcb_struct * polyder, long len, long prec)
{
    fmpcb_t t;
    fmpr_t u, v;

    fmpcb_init(t);
    fmpr_init(u);
    fmpr_init(v);

    fmpcb_set(r, m);
    fmpr_zero(fmprb_radref(fmpcb_realref(r)));
    fmpr_zero(fmprb_radref(fmpcb_imagref(r)));

    _fmpcb_poly_evaluate(t, poly, len, r, prec);
    fmpcb_get_abs_ubound_fmpr(u, t, FMPRB_RAD_PREC);

    /* it could happen that we have an exact root, in which case
       we should avoid dividing by the derivative */
    if (!fmpr_is_zero(u))
    {
        _fmpcb_poly_evaluate(t, polyder, len - 1, r, prec);
        fmpcb_inv(t, t, FMPRB_RAD_PREC);
        fmpcb_get_abs_ubound_fmpr(v, t, FMPRB_RAD_PREC);

        fmpr_mul(u, u, v, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul_ui(u, u, len - 1, FMPRB_RAD_PREC, FMPR_RND_UP);
    }

    fmpr_set(fmprb_radref(fmpcb_realref(r)), u);
    fmpr_set(fmprb_radref(fmpcb_imagref(r)), u);

    fmpr_clear(u);
    fmpr_clear(v);
    fmpcb_clear(t);
}

