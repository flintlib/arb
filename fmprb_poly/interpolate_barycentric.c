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

void
_fmprb_poly_interpolate_barycentric(fmprb_ptr poly,
    fmprb_srcptr xs, fmprb_srcptr ys, long n, long prec)
{
    fmprb_ptr P, Q, w;
    fmprb_t t;
    long i, j;

    if (n == 1)
    {
        fmprb_set(poly, ys);
        return;
    }

    P = _fmprb_vec_init(n + 1);
    Q = _fmprb_vec_init(n);
    w = _fmprb_vec_init(n);
    fmprb_init(t);

    _fmprb_poly_product_roots(P, xs, n, prec);

    for (i = 0; i < n; i++)
    {
        fmprb_one(w + i);

        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                fmprb_sub(t, xs + i, xs + j, prec);
                fmprb_mul(w + i, w + i, t, prec);
            }
        }

        fmprb_ui_div(w + i, 1UL, w + i, prec);
    }

    _fmprb_vec_zero(poly, n);

    for (i = 0; i < n; i++)
    {
        _fmprb_poly_div_root(Q, t, P, n + 1, xs + i, prec);
        fmprb_mul(t, w + i, ys + i, prec);
        _fmprb_vec_scalar_addmul(poly, Q, n, t, prec);
    }

    _fmprb_vec_clear(P, n + 1);
    _fmprb_vec_clear(Q, n);
    _fmprb_vec_clear(w, n);
    fmprb_clear(t);
}

void
fmprb_poly_interpolate_barycentric(fmprb_poly_t poly,
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
        _fmprb_poly_interpolate_barycentric(poly->coeffs,
            xs, ys, n, prec);
        _fmprb_poly_normalise(poly);
    }
}
