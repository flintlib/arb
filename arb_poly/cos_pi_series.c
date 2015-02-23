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
_arb_poly_cos_pi_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        arb_cos_pi(g, h, prec);
        _arb_vec_zero(g + 1, n - 1);
    }
    else if (n == 2)
    {
        arb_t t;
        arb_init(t);
        arb_sin_cos_pi(t, g, h, prec);
        arb_neg(t, t);
        arb_mul(g + 1, h + 1, t, prec);  /* safe since hlen >= 2 */
        arb_const_pi(t, prec);
        arb_mul(g + 1, g + 1, t, prec);
        arb_clear(t);
    }
    else
    {
        arb_ptr t = _arb_vec_init(n);
        _arb_poly_sin_cos_pi_series(t, g, h, hlen, n, prec);
        _arb_vec_clear(t, n);
    }
}

void
arb_poly_cos_pi_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    if (hlen == 0)
    {
        arb_poly_one(g);
        return;
    }

    if (hlen == 1)
        n = 1;

    arb_poly_fit_length(g, n);
    _arb_poly_cos_pi_series(g->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}

