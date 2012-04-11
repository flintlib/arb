/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

void
_arb_poly_set_coeff_same_exp(arb_poly_t z, const arb_t c)
{
    long e1, e2;
    fmpz_t rad;

    e1 = z->exp;
    e2 = c->exp;

    fmpz_init(rad);

    if (e1 > e2)
    {
        fmpz_tdiv_q_2exp(z->coeffs, arb_midref(c), e1 - e2);
        fmpz_cdiv_q_2exp(rad, arb_radref(c), e1 - e2);
        fmpz_add_ui(rad, rad, 1);
    }
    else
    {
        fmpz_mul_2exp(z->coeffs, arb_midref(c), e2 - e1);
        fmpz_mul_2exp(rad, arb_radref(c), e2 - e1);
    }

    if (fmpz_cmp(arb_poly_radref(z), rad) < 0)
        fmpz_set(arb_poly_radref(z), rad);

    fmpz_clear(rad);
}

void
arb_poly_log_series(arb_poly_t z, const arb_poly_t x, long n)
{
    arb_poly_t t, u;
    arb_t c;

    arb_poly_init(t, arb_poly_prec(z));
    arb_poly_init(u, arb_poly_prec(z));
    arb_init(c, arb_poly_prec(z));

    /* first term */
    fmpz_set(arb_midref(c), x->coeffs);
    fmpz_set(arb_radref(c), arb_poly_radref(x));
    fmpz_set(arb_expref(c), arb_poly_expref(x));
    arb_log(c, c);

    arb_poly_inv_series(t, x, n);
    arb_poly_derivative(u, x);
    arb_poly_mullow(t, t, u, n - 1);
    arb_poly_integral(z, t);

    _arb_poly_set_coeff_same_exp(z, c);

    arb_clear(c);
    arb_poly_clear(t);
    arb_poly_clear(u);
}
