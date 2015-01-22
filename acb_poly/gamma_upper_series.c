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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
_acb_poly_gamma_upper_series(acb_ptr g, const acb_t s, acb_srcptr h, long hlen, long n, long prec)
{
    acb_t c;
    acb_init(c);

    acb_hypgeom_gamma_upper(c, s, h, 0, prec);

    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        _acb_vec_zero(g + 1, n - 1);
    }
    else
    {
        acb_ptr t, u, v;

        t = _acb_vec_init(n);
        u = _acb_vec_init(n);
        v = _acb_vec_init(n);

        /* Gamma(s, h(x)) = -integral(h'(x) h(x)^(s-1) exp(-h(x)) */
        acb_sub_ui(u, s, 1, prec);
        _acb_poly_pow_acb_series(t, h, hlen, u, n, prec);
        _acb_poly_derivative(u, h, hlen, prec);
        _acb_poly_mullow(v, t, n, u, hlen - 1, n, prec);
        _acb_vec_neg(t, h, hlen);
        _acb_poly_exp_series(t, t, hlen, n, prec);
        _acb_poly_mullow(g, v, n, t, n, n, prec);
        _acb_poly_integral(g, g, n, prec);
        _acb_vec_neg(g, g, n);

        _acb_vec_clear(t, n);
        _acb_vec_clear(u, n);
        _acb_vec_clear(v, n);
    }

    acb_swap(g, c);
    acb_clear(c);
}

void
acb_poly_gamma_upper_series(acb_poly_t g, const acb_t s, const acb_poly_t h, long n, long prec)
{
    long hlen = h->length;

    if (n == 0)
    {
        acb_poly_zero(g);
        return;
    }

    acb_poly_fit_length(g, n);

    if (hlen == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_poly_gamma_upper_series(g->coeffs, s, t, 1, n, prec);
        acb_clear(t);
    }
    else
    {
        _acb_poly_gamma_upper_series(g->coeffs, s, h->coeffs, hlen, n, prec);
    }

    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}

