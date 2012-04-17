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
arb_poly_exp_series(arb_poly_t z, const arb_poly_t x, long n, int exact_const)
{
    long a[FLINT_BITS];
    long i;
    arb_poly_t t, v, c;
    arb_t r;

    if (n == 0 || x->length == 0)
    {
        arb_poly_zero(z);
        return;
    }

    a[i = 0] = n;
    while (n >= 2)
        a[++i] = (n = (n + 1) / 2);

    arb_poly_init(t, arb_poly_prec(z));

    /* remove constant term */
    arb_poly_init(v, arb_poly_prec(x));
    arb_poly_set(v, x);
    fmpz_zero(arb_poly_coeffs(v));

    arb_poly_init(c, arb_poly_prec(z));

    if (!exact_const)
    {
        /* first coefficient */
        arb_init(r, arb_poly_prec(z));
        fmpz_set(arb_midref(r), arb_poly_coeffs(x));
        fmpz_set(arb_radref(r), arb_poly_radref(x));
        fmpz_set(arb_expref(r), arb_poly_expref(x));

        arb_exp(r, r);
        _arb_poly_fit_length(c, 1);
        fmpz_set(arb_poly_coeffs(c), arb_midref(r));
        fmpz_set(arb_poly_radref(c), arb_radref(r));
        fmpz_set(arb_poly_expref(c), arb_expref(r));
        c->length = 1;
        arb_clear(r);
    }

    arb_poly_set_si(z, 1);

    for (i--; i >= 0; i--)
    {
        n = a[i];

        arb_poly_log_series(t, z, n, 1);
        arb_poly_neg(t, t);
        arb_poly_add(t, t, v);
        arb_poly_mullow(t, z, t, n);
        arb_poly_add(z, z, t);
    }

    if (!exact_const)
        arb_poly_mul(z, z, c);

    arb_poly_clear(t);
    arb_poly_clear(v);
    arb_poly_clear(c);
}
