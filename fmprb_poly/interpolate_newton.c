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

#include "fmprb_poly.h"

static void
_interpolate_newton(fmprb_ptr ys, fmprb_srcptr xs, long n, long prec)
{
    fmprb_t p, q, t;
    long i, j;

    fmprb_init(p);
    fmprb_init(q);
    fmprb_init(t);

    for (i = 1; i < n; i++)
    {
        fmprb_set(t, ys + i - 1);

        for (j = i; j < n; j++)
        {
            fmprb_sub(p, ys + j, t, prec);
            fmprb_sub(q, xs + j, xs + j - i, prec);
            fmprb_set(t, ys + j);
            fmprb_div(ys + j, p, q, prec);
        }
    }

    fmprb_clear(p);
    fmprb_clear(q);
    fmprb_clear(t);
}

static void
_newton_to_monomial(fmprb_ptr ys, fmprb_srcptr xs, long n, long prec)
{
    fmprb_t t, u;
    long i, j;

    fmprb_init(t);
    fmprb_init(u);

    for (i = n - 2; i >= 0; i--)
    {
        fmprb_set(t, ys + i);
        fmprb_set(ys + i, ys + i + 1);

        for (j = i + 1; j < n - 1; j++)
        {
            fmprb_mul(u, ys + j, xs + i, prec);
            fmprb_sub(ys + j, ys + j + 1, u, prec);
        }

        fmprb_mul(u, ys + n - 1, xs + i, prec);
        fmprb_sub(ys + n - 1, t, u, prec);
    }

    _fmprb_poly_reverse(ys, ys, n, n);

    fmprb_clear(t);
    fmprb_clear(u);
}

void
_fmprb_poly_interpolate_newton(fmprb_ptr poly, fmprb_srcptr xs,
    fmprb_srcptr ys, long n, long prec)
{
    if (n == 1)
    {
        fmprb_set(poly, ys);
    }
    else
    {
        _fmprb_vec_set(poly, ys, n);
        _interpolate_newton(poly, xs, n, prec);
        while (n > 0 && fmprb_is_zero(poly + n - 1)) n--;
        _newton_to_monomial(poly, xs, n, prec);
    }
}

void
fmprb_poly_interpolate_newton(fmprb_poly_t poly,
    fmprb_srcptr xs, fmprb_srcptr ys, long n, long prec)
{
    if (n == 0)
    {
        fmprb_poly_zero(poly);
    }
    else
    {
        fmprb_poly_fit_length(poly, n);
        _fmprb_poly_set_length(poly, n);
        _fmprb_poly_interpolate_newton(poly->coeffs,
            xs, ys, n, prec);
        _fmprb_poly_normalise(poly);
    }
}
