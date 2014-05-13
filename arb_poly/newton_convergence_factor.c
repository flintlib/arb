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

#include "arb_poly.h"

void
_arb_poly_newton_convergence_factor(arf_t convergence_factor,
    arb_srcptr poly, long len,
    const arb_t convergence_interval, long prec)
{
    arb_ptr deriv;
    arb_t t, u;

    arb_init(t);
    arb_init(u);
    deriv = _arb_vec_init(len - 1);

    _arb_poly_derivative(deriv, poly, len, prec);
    _arb_poly_evaluate(t, deriv, len - 1, convergence_interval, prec);

    _arb_poly_derivative(deriv, deriv, len - 1, prec);
    _arb_poly_evaluate(u, deriv, len - 2, convergence_interval, prec);

    arb_div(t, u, t, prec);
    arb_mul_2exp_si(t, t, -1);

    arb_get_abs_ubound_arf(convergence_factor, t, prec);

    _arb_vec_clear(deriv, len - 1);
    arb_clear(t);
    arb_clear(u);
}

