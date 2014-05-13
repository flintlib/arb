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

#include "arb_poly.h"

void
_arb_poly_interpolation_weights(arb_ptr w,
    arb_ptr * tree, long len, long prec)
{
    arb_ptr tmp;
    long i, n, height;

    if (len == 0)
        return;

    if (len == 1)
    {
        arb_one(w);
        return;
    }

    tmp = _arb_vec_init(len + 1);
    height = FLINT_CLOG2(len);
    n = 1L << (height - 1);

    _arb_poly_mul_monic(tmp, tree[height-1], n + 1,
                        tree[height-1] + (n + 1), (len - n + 1), prec);

    _arb_poly_derivative(tmp, tmp, len + 1, prec);
    _arb_poly_evaluate_vec_fast_precomp(w, tmp, len, tree, len, prec);

    for (i = 0; i < len; i++)
        arb_inv(w + i, w + i, prec);

    _arb_vec_clear(tmp, len + 1);
}

void
_arb_poly_interpolate_fast_precomp(arb_ptr poly,
    arb_srcptr ys, arb_ptr * tree, arb_srcptr weights,
    long len, long prec)
{
    arb_ptr t, u, pa, pb;
    long i, pow, left;

    if (len == 0)
        return;

    t = _arb_vec_init(len);
    u = _arb_vec_init(len);

    for (i = 0; i < len; i++)
        arb_mul(poly + i, weights + i, ys + i, prec);

    for (i = 0; i < FLINT_CLOG2(len); i++)
    {
        pow = (1L << i);
        pa = tree[i];
        pb = poly;
        left = len;

        while (left >= 2 * pow)
        {
            _arb_poly_mul(t, pa, pow + 1, pb + pow, pow, prec);
            _arb_poly_mul(u, pa + pow + 1, pow + 1, pb, pow, prec);
            _arb_vec_add(pb, t, u, 2 * pow, prec);

            left -= 2 * pow;
            pa += 2 * pow + 2;
            pb += 2 * pow;
        }

        if (left > pow)
        {
            _arb_poly_mul(t, pa, pow + 1, pb + pow, left - pow, prec);
            _arb_poly_mul(u, pb, pow, pa + pow + 1, left - pow + 1, prec);
            _arb_vec_add(pb, t, u, left, prec);
        }
    }

    _arb_vec_clear(t, len);
    _arb_vec_clear(u, len);
}

void
_arb_poly_interpolate_fast(arb_ptr poly,
    arb_srcptr xs, arb_srcptr ys, long len, long prec)
{
    arb_ptr * tree;
    arb_ptr w;

    tree = _arb_poly_tree_alloc(len);
    _arb_poly_tree_build(tree, xs, len, prec);

    w = _arb_vec_init(len);
    _arb_poly_interpolation_weights(w, tree, len, prec);

    _arb_poly_interpolate_fast_precomp(poly, ys, tree, w, len, prec);

    _arb_vec_clear(w, len);
    _arb_poly_tree_free(tree, len);
}

void
arb_poly_interpolate_fast(arb_poly_t poly,
        arb_srcptr xs, arb_srcptr ys, long n, long prec)
{
    if (n == 0)
    {
        arb_poly_zero(poly);
    }
    else
    {
        arb_poly_fit_length(poly, n);
        _arb_poly_set_length(poly, n);
        _arb_poly_interpolate_fast(poly->coeffs, xs, ys, n, prec);
        _arb_poly_normalise(poly);
    }
}
