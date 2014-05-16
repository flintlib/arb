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

acb_ptr * _acb_poly_tree_alloc(long len)
{
    acb_ptr * tree = NULL;

    if (len)
    {
        long i, height = FLINT_CLOG2(len);

        tree = flint_malloc(sizeof(acb_ptr) * (height + 1));
        for (i = 0; i <= height; i++)
            tree[i] = _acb_vec_init(len + (len >> i) + 1);
    }

    return tree;
}

void _acb_poly_tree_free(acb_ptr * tree, long len)
{
    if (len)
    {
        long i, height = FLINT_CLOG2(len);

        for (i = 0; i <= height; i++)
            _acb_vec_clear(tree[i], len + (len >> i) + 1);

        flint_free(tree);
    }
}

void
_acb_poly_tree_build(acb_ptr * tree, acb_srcptr roots, long len, long prec)
{
    long height, pow, left, i;
    acb_ptr pa, pb;
    acb_srcptr a, b;

    if (len == 0)
        return;

    height = FLINT_CLOG2(len);

    /* zeroth level, (x-a) */
    for (i = 0; i < len; i++)
    {
        acb_one(tree[0] + (2 * i + 1));
        acb_neg(tree[0] + (2 * i), roots + i);
    }

    /* first level, (x-a)(x-b) = x^2 + (-a-b)*x + a*b */
    if (height > 1)
    {
        pa = tree[1];

        for (i = 0; i < len / 2; i++)
        {
            a = (acb_srcptr) (roots + (2 * i));
            b = (acb_srcptr) (roots + (2 * i + 1));

            acb_mul(pa + (3 * i), a, b, prec);
            acb_add(pa + (3 * i + 1), a, b, prec);
            acb_neg(pa + (3 * i + 1), pa + (3 * i + 1));
            acb_one(pa + (3 * i + 2));
        }

        if (len & 1)
        {
            acb_neg(pa + (3 * (len / 2)), roots + len - 1);
            acb_one(pa + (3 * (len / 2) + 1));
        }
    }

    for (i = 1; i < height - 1; i++)
    {
        left = len;
        pow = 1L << i;
        pa = tree[i];
        pb = tree[i + 1];

        while (left >= 2 * pow)
        {
            _acb_poly_mul_monic(pb, pa, pow + 1, pa + pow + 1, pow + 1, prec);
            left -= 2 * pow;
            pa += 2 * pow + 2;
            pb += 2 * pow + 1;
        }

        if (left > pow)
        {
            _acb_poly_mul_monic(pb, pa, pow + 1, pa + pow + 1, left - pow + 1, prec);
        }
        else if (left > 0)
            _acb_vec_set(pb, pa, left + 1);
    }
}
