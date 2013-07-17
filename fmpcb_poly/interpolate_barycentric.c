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

void
_fmpcb_poly_interpolate_barycentric(fmpcb_ptr poly,
    fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec)
{
    fmpcb_ptr P, Q, w;
    fmpcb_t t;
    long i, j;

    if (n == 1)
    {
        fmpcb_set(poly, ys);
        return;
    }

    P = _fmpcb_vec_init(n + 1);
    Q = _fmpcb_vec_init(n);
    w = _fmpcb_vec_init(n);
    fmpcb_init(t);

    _fmpcb_poly_product_roots(P, xs, n, prec);

    for (i = 0; i < n; i++)
    {
        fmpcb_one(w + i);

        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                fmpcb_sub(t, xs + i, xs + j, prec);
                fmpcb_mul(w + i, w + i, t, prec);
            }
        }

        fmpcb_inv(w + i, w + i, prec);
    }

    _fmpcb_vec_zero(poly, n);

    for (i = 0; i < n; i++)
    {
        _fmpcb_poly_div_root(Q, t, P, n + 1, xs + i, prec);
        fmpcb_mul(t, w + i, ys + i, prec);
        _fmpcb_vec_scalar_addmul(poly, Q, n, t, prec);
    }

    _fmpcb_vec_clear(P, n + 1);
    _fmpcb_vec_clear(Q, n);
    _fmpcb_vec_clear(w, n);
    fmpcb_clear(t);
}

void
fmpcb_poly_interpolate_barycentric(fmpcb_poly_t poly,
    fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec)
{
    if (n == 0)
    {
        fmpcb_poly_zero(poly);
    }
    else
    {
        fmpcb_poly_fit_length(poly, n);
        _fmpcb_poly_set_length(poly, n);
        _fmpcb_poly_interpolate_barycentric(poly->coeffs,
            xs, ys, n, prec);
        _fmpcb_poly_normalise(poly);
    }
}
