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
fmpcb_sqrt(fmpcb_t y, const fmpcb_t x, long prec)
{
    fmprb_t r, t, u;
    long wp;

#define a fmpcb_realref(x)
#define b fmpcb_imagref(x)
#define c fmpcb_realref(y)
#define d fmpcb_imagref(y)

    if (fmprb_is_zero(b))
    {
        if (fmprb_is_nonnegative(a))
        {
            fmprb_sqrt(c, a, prec);
            fmprb_zero(d);
            return;
        }
        else if (fmprb_is_nonpositive(a))
        {
            fmprb_neg(d, a);
            fmprb_sqrt(d, d, prec);
            fmprb_zero(c);
            return;
        }
    }

    if (fmprb_is_zero(a))
    {
        if (fmprb_is_nonnegative(b))
        {
            fmprb_mul_2exp_si(c, b, -1);
            fmprb_sqrt(c, c, prec);
            fmprb_set(d, c);
            return;
        }
        else if (fmprb_is_nonpositive(b))
        {
            fmprb_mul_2exp_si(c, b, -1);
            fmprb_neg(c, c);
            fmprb_sqrt(c, c, prec);
            fmprb_neg(d, c);
            return;
        }
    }

    /* sqrt(a+bi) = sqrt((r+a)/2) + b/sqrt(2*(r+a))*i, r = |a+bi| */

    wp = prec + 4;

    fmprb_init(r);
    fmprb_init(t);
    fmprb_init(u);

    fmpcb_abs(r, x, wp);
    fmprb_add(t, r, a, wp);

    fmprb_mul_2exp_si(u, t, 1);
    fmprb_sqrt(u, u, wp);
    fmprb_div(d, b, u, prec);

    fmprb_set_round(c, u, prec);
    fmprb_mul_2exp_si(c, c, -1);

    fmprb_clear(r);
    fmprb_clear(t);
    fmprb_clear(u);

#undef a
#undef b
#undef c
#undef d
}

