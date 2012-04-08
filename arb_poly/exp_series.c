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
arb_poly_exp_series(arb_poly_t z, const arb_poly_t x, long n)
{
    long a[FLINT_BITS];
    long i, m;
    arb_poly_t t;
    long xlen;

    xlen = x->length;

    if (z == x)
    {
        arb_poly_t t;
        arb_poly_init(t, arb_poly_prec(x));
        arb_poly_set(t, x);
        arb_poly_exp_series(z, t, n);
        arb_poly_clear(t);
        return;
    }

    a[i = 0] = n;
    while (n >= 2)
        a[++i] = (n = (n + 1) / 2);

    arb_poly_init(t, arb_poly_prec(z));

    /* first coefficient (assuming input zero, todo: general case) */
    arb_poly_set_si(z, 1);
    fmpz_mul_2exp(arb_poly_coeffs(z), arb_poly_coeffs(z), arb_poly_prec(z));
    fmpz_set_si(arb_poly_expref(z), -arb_poly_prec(z));
    fmpz_zero(arb_poly_radref(z));

    /* Newton iteration */
    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];

        arb_poly_log_series(t, z, n);
        arb_poly_neg(t, t);
        arb_poly_add(t, t, x);
        arb_poly_mullow(t, z, t, n);
        arb_poly_add(z, z, t);
    }

    arb_poly_clear(t);
}
