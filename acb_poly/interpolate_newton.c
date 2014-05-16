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

#include "acb_poly.h"

static void
_interpolate_newton(acb_ptr ys, acb_srcptr xs, long n, long prec)
{
    acb_t p, q, t;
    long i, j;

    acb_init(p);
    acb_init(q);
    acb_init(t);

    for (i = 1; i < n; i++)
    {
        acb_set(t, ys + i - 1);

        for (j = i; j < n; j++)
        {
            acb_sub(p, ys + j, t, prec);
            acb_sub(q, xs + j, xs + j - i, prec);
            acb_set(t, ys + j);
            acb_div(ys + j, p, q, prec);
        }
    }

    acb_clear(p);
    acb_clear(q);
    acb_clear(t);
}

static void
_newton_to_monomial(acb_ptr ys, acb_srcptr xs, long n, long prec)
{
    acb_t t, u;
    long i, j;

    acb_init(t);
    acb_init(u);

    for (i = n - 2; i >= 0; i--)
    {
        acb_set(t, ys + i);
        acb_set(ys + i, ys + i + 1);

        for (j = i + 1; j < n - 1; j++)
        {
            acb_mul(u, ys + j, xs + i, prec);
            acb_sub(ys + j, ys + j + 1, u, prec);
        }

        acb_mul(u, ys + n - 1, xs + i, prec);
        acb_sub(ys + n - 1, t, u, prec);
    }

    _acb_poly_reverse(ys, ys, n, n);

    acb_clear(t);
    acb_clear(u);
}

void
_acb_poly_interpolate_newton(acb_ptr poly, acb_srcptr xs,
    acb_srcptr ys, long n, long prec)
{
    if (n == 1)
    {
        acb_set(poly, ys);
    }
    else
    {
        _acb_vec_set(poly, ys, n);
        _interpolate_newton(poly, xs, n, prec);
        while (n > 0 && acb_is_zero(poly + n - 1)) n--;
        _newton_to_monomial(poly, xs, n, prec);
    }
}

void
acb_poly_interpolate_newton(acb_poly_t poly,
    acb_srcptr xs, acb_srcptr ys, long n, long prec)
{
    if (n == 0)
    {
        acb_poly_zero(poly);
    }
    else
    {
        acb_poly_fit_length(poly, n);
        _acb_poly_set_length(poly, n);
        _acb_poly_interpolate_newton(poly->coeffs,
            xs, ys, n, prec);
        _acb_poly_normalise(poly);
    }
}
