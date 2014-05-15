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
acb_sqrt(acb_t y, const acb_t x, long prec)
{
    arb_t r, t, u;
    long wp;

#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(y)
#define d acb_imagref(y)

    if (arb_is_zero(b))
    {
        if (arb_is_nonnegative(a))
        {
            arb_sqrt(c, a, prec);
            arb_zero(d);
            return;
        }
        else if (arb_is_nonpositive(a))
        {
            arb_neg(d, a);
            arb_sqrt(d, d, prec);
            arb_zero(c);
            return;
        }
    }

    if (arb_is_zero(a))
    {
        if (arb_is_nonnegative(b))
        {
            arb_mul_2exp_si(c, b, -1);
            arb_sqrt(c, c, prec);
            arb_set(d, c);
            return;
        }
        else if (arb_is_nonpositive(b))
        {
            arb_mul_2exp_si(c, b, -1);
            arb_neg(c, c);
            arb_sqrt(c, c, prec);
            arb_neg(d, c);
            return;
        }
    }

    /* sqrt(a+bi) = sqrt((r+a)/2) + b/sqrt(2*(r+a))*i, r = |a+bi| */

    wp = prec + 4;

    arb_init(r);
    arb_init(t);
    arb_init(u);

    acb_abs(r, x, wp);
    arb_add(t, r, a, wp);

    arb_mul_2exp_si(u, t, 1);
    arb_sqrt(u, u, wp);
    arb_div(d, b, u, prec);

    arb_set_round(c, u, prec);
    arb_mul_2exp_si(c, c, -1);

    arb_clear(r);
    arb_clear(t);
    arb_clear(u);

#undef a
#undef b
#undef c
#undef d
}

