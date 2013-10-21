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

#include "fmpcb.h"

void
fmpcb_rsqrt(fmpcb_t y, const fmpcb_t x, long prec)
{
    fmprb_t r, t, u, v;
    long wp;

#define a fmpcb_realref(x)
#define b fmpcb_imagref(x)
#define c fmpcb_realref(y)
#define d fmpcb_imagref(y)

    if (fmprb_is_zero(b))
    {
        if (fmprb_is_nonnegative(a))
        {
            fmprb_rsqrt(c, a, prec);
            fmprb_zero(d);
            return;
        }
        else if (fmprb_is_nonpositive(a))
        {
            fmprb_neg(d, a);
            fmprb_rsqrt(d, d, prec);
            fmprb_neg(d, d);
            fmprb_zero(c);
            return;
        }
    }

    if (fmprb_is_zero(a))
    {
        if (fmprb_is_nonnegative(b))
        {
            fmprb_mul_2exp_si(c, b, 1);
            fmprb_rsqrt(c, c, prec);
            fmprb_neg(d, c);
            return;
        }
        else if (fmprb_is_nonpositive(b))
        {
            fmprb_mul_2exp_si(c, b, 1);
            fmprb_neg(c, c);
            fmprb_rsqrt(c, c, prec);
            fmprb_set(d, c);
            return;
        }
    }

    /* based on the identity sqrt(z) = sqrt(r) (z+r) / |z+r|: */
    /* 1/sqrt(a+bi) = (1/v)((a+r) - b*i), r = |a+bi|, v = sqrt(r*(b^2+(a+r)^2)) */

    wp = prec + 6;

    fmprb_init(r);
    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    /* u = b^2, r = |a+bi| */
    fmprb_mul(t, a, a, wp);
    fmprb_mul(u, b, b, wp);
    fmprb_add(r, t, u, wp);
    fmprb_sqrtpos(r, r, wp);

    /* t = a+r, v = r*(b^2+(a+r)^2) */
    fmprb_add(t, r, a, wp);
    fmprb_mul(v, t, t, wp);
    fmprb_add(v, v, u, wp);
    fmprb_mul(v, v, r, wp);

    /* v = 1/sqrt(v) */
    fmprb_rsqrt(v, v, wp);

    fmprb_mul(c, t, v, prec);
    fmprb_mul(d, b, v, prec);
    fmprb_neg(d, d);

    fmprb_clear(r);
    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);

#undef a
#undef b
#undef c
#undef d
}

