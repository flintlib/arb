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

#include "acb.h"

void
acb_rsqrt(acb_t y, const acb_t x, long prec)
{
    arb_t r, t, u, v;
    long wp;

#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(y)
#define d acb_imagref(y)

    if (arb_is_zero(b))
    {
        if (arb_is_nonnegative(a))
        {
            arb_rsqrt(c, a, prec);
            arb_zero(d);
            return;
        }
        else if (arb_is_nonpositive(a))
        {
            arb_neg(d, a);
            arb_rsqrt(d, d, prec);
            arb_neg(d, d);
            arb_zero(c);
            return;
        }
    }

    if (arb_is_zero(a))
    {
        if (arb_is_nonnegative(b))
        {
            arb_mul_2exp_si(c, b, 1);
            arb_rsqrt(c, c, prec);
            arb_neg(d, c);
            return;
        }
        else if (arb_is_nonpositive(b))
        {
            arb_mul_2exp_si(c, b, 1);
            arb_neg(c, c);
            arb_rsqrt(c, c, prec);
            arb_set(d, c);
            return;
        }
    }

    /* based on the identity sqrt(z) = sqrt(r) (z+r) / |z+r|: */
    /* 1/sqrt(a+bi) = (1/v)((a+r) - b*i), r = |a+bi|, v = sqrt(r*(b^2+(a+r)^2)) */

    wp = prec + 6;

    arb_init(r);
    arb_init(t);
    arb_init(u);
    arb_init(v);

    /* u = b^2, r = |a+bi| */
    arb_mul(t, a, a, wp);
    arb_mul(u, b, b, wp);
    arb_add(r, t, u, wp);
    arb_sqrtpos(r, r, wp);

    /* t = a+r, v = r*(b^2+(a+r)^2) */
    arb_add(t, r, a, wp);
    arb_mul(v, t, t, wp);
    arb_add(v, v, u, wp);
    arb_mul(v, v, r, wp);

    /* v = 1/sqrt(v) */
    arb_rsqrt(v, v, wp);

    arb_mul(c, t, v, prec);
    arb_mul(d, b, v, prec);
    arb_neg(d, d);

    arb_clear(r);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);

#undef a
#undef b
#undef c
#undef d
}

