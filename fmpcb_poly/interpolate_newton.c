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

#include "fmpcb_poly.h"

static void
_interpolate_newton(fmpcb_struct * ys, const fmpcb_struct * xs, long n, long prec)
{
    fmpcb_t p, q, t;
    long i, j;

    fmpcb_init(p);
    fmpcb_init(q);
    fmpcb_init(t);

    for (i = 1; i < n; i++)
    {
        fmpcb_set(t, ys + i - 1);

        for (j = i; j < n; j++)
        {
            fmpcb_sub(p, ys + j, t, prec);
            fmpcb_sub(q, xs + j, xs + j - i, prec);
            fmpcb_set(t, ys + j);
            fmpcb_div(ys + j, p, q, prec);
        }
    }

    fmpcb_clear(p);
    fmpcb_clear(q);
    fmpcb_clear(t);
}

static void
_newton_to_monomial(fmpcb_struct * ys, const fmpcb_struct * xs, long n, long prec)
{
    fmpcb_t t, u;
    long i, j;

    fmpcb_init(t);
    fmpcb_init(u);

    for (i = n - 2; i >= 0; i--)
    {
        fmpcb_set(t, ys + i);
        fmpcb_set(ys + i, ys + i + 1);

        for (j = i + 1; j < n - 1; j++)
        {
            fmpcb_mul(u, ys + j, xs + i, prec);
            fmpcb_sub(ys + j, ys + j + 1, u, prec);
        }

        fmpcb_mul(u, ys + n - 1, xs + i, prec);
        fmpcb_sub(ys + n - 1, t, u, prec);
    }

    _fmpcb_poly_reverse(ys, ys, n, n);

    fmpcb_clear(t);
    fmpcb_clear(u);
}

void
_fmpcb_poly_interpolate_newton(fmpcb_struct * poly, const fmpcb_struct * xs,
    const fmpcb_struct * ys, long n, long prec)
{
    if (n == 1)
    {
        fmpcb_set(poly, ys);
    }
    else
    {
        _fmpcb_vec_set(poly, ys, n);
        _interpolate_newton(poly, xs, n, prec);
        while (n > 0 && fmpcb_is_zero(poly + n - 1)) n--;
        _newton_to_monomial(poly, xs, n, prec);
    }
}

void
fmpcb_poly_interpolate_newton(fmpcb_poly_t poly,
    const fmpcb_struct * xs, const fmpcb_struct * ys, long n, long prec)
{
    if (n == 0)
    {
        fmpcb_poly_zero(poly);
    }
    else
    {
        fmpcb_poly_fit_length(poly, n);
        _fmpcb_poly_set_length(poly, n);
        _fmpcb_poly_interpolate_newton(poly->coeffs,
            xs, ys, n, prec);
        _fmpcb_poly_normalise(poly);
    }
}
