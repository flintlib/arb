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

void
_fmprb_poly_newton_convergence_factor(fmpr_t convergence_factor,
    fmprb_srcptr poly, long len,
    const fmprb_t convergence_interval, long prec)
{
    fmprb_ptr deriv;
    fmprb_t t, u;

    fmprb_init(t);
    fmprb_init(u);
    deriv = _fmprb_vec_init(len - 1);

    _fmprb_poly_derivative(deriv, poly, len, prec);
    _fmprb_poly_evaluate(t, deriv, len - 1, convergence_interval, prec);

    _fmprb_poly_derivative(deriv, deriv, len - 1, prec);
    _fmprb_poly_evaluate(u, deriv, len - 2, convergence_interval, prec);

    fmprb_div(t, u, t, prec);
    fmprb_mul_2exp_si(t, t, -1);

    fmprb_get_abs_ubound_fmpr(convergence_factor, t, prec);

    _fmprb_vec_clear(deriv, len - 1);
    fmprb_clear(t);
    fmprb_clear(u);
}

