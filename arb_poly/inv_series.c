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

/* todo: fix arb_poly_add so this isn't needed */
static void
_arb_poly_add_two(arb_poly_t t)
{
    if (fmpz_sgn(arb_poly_expref(t)) <= 0)
    {
        fmpz_t two;
        fmpz_init(two);
        fmpz_set_ui(two, 2UL);
        fmpz_mul_2exp(two, two, -fmpz_get_si(arb_poly_expref(t)));
        fmpz_add(t->coeffs, t->coeffs, two);
        fmpz_clear(two);
    }
    else if (*arb_poly_expref(t) == 1)
    {
        fmpz_add_ui(t->coeffs, t->coeffs, 1);
    }
    else
    {
        fmpz_add_ui(arb_poly_radref(t), arb_poly_radref(t), 1);
    }
}

void
arb_poly_inv_series(arb_poly_t z, const arb_poly_t x, long n)
{
    long a[FLINT_BITS];
    long i;
    arb_poly_t t, u;
    arb_t one, r;
    long xlen;

    xlen = x->length;

    if (xlen == 0)
    {
        printf("arb_poly_inv_series: division by zero\n");
        abort();
    }

    if (z == x)
    {
        arb_poly_t t;
        arb_poly_init(t, arb_poly_prec(x));
        arb_poly_set(t, x);
        arb_poly_inv_series(z, t, n);
        arb_poly_clear(t);
        return;
    }

    a[i = 0] = n;
    while (n >= 2)
        a[++i] = (n = (n + 1) / 2);

    arb_poly_init(t, arb_poly_prec(z));
    arb_poly_init(u, arb_poly_prec(z));

    /* first coefficient */
    arb_init(r, arb_poly_prec(z));
    arb_init(one, arb_poly_prec(z));
    arb_set_si(one, 1);

    fmpz_set(arb_midref(r), arb_poly_coeffs(x));
    fmpz_set(arb_radref(r), arb_poly_radref(x));
    fmpz_set(arb_expref(r), arb_poly_expref(x));

    arb_div(r, one, r);

    _arb_poly_fit_length(z, 1);
    fmpz_set(arb_poly_coeffs(z), arb_midref(r));
    fmpz_set(arb_poly_radref(z), arb_radref(r));
    fmpz_set(arb_poly_expref(z), arb_expref(r));
    z->length = 1;

    arb_clear(r);
    arb_clear(one);
    /* done first coefficient */

    /* Newton iteration */
    for (i--; i >= 0; i--)
    {
        n = a[i];

        arb_poly_mullow(t, x, z, n);
        arb_poly_neg(t, t);

        _arb_poly_add_two(t);

        arb_poly_mullow(u, t, z, n);
        arb_poly_set(z, u);
    }

    arb_poly_clear(t);
    arb_poly_clear(u);
}
