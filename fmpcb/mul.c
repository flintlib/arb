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

#include "fmpcb.h"

void
fmpcb_mul(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)
{
#define a fmpcb_realref(x)
#define b fmpcb_imagref(x)
#define c fmpcb_realref(y)
#define d fmpcb_imagref(y)
#define e fmpcb_realref(z)
#define f fmpcb_imagref(z)

    if (fmprb_is_zero(b))
    {
        fmprb_mul(f, d, a, prec);
        fmprb_mul(e, c, a, prec);
    }
    else if (fmprb_is_zero(d))
    {
        fmprb_mul(f, b, c, prec);
        fmprb_mul(e, a, c, prec);
    }
    else if (fmprb_is_zero(a))
    {
        fmprb_mul(e, c, b, prec);
        fmprb_mul(f, d, b, prec);
        fmpcb_mul_onei(z, z);
    }
    else if (fmprb_is_zero(c))
    {
        fmprb_mul(e, a, d, prec);
        fmprb_mul(f, b, d, prec);
        fmpcb_mul_onei(z, z);
    }
    /* squaring = a^2-b^2, 2ab */
    else if (x == y)
    {
        /* aliasing */
        if (z == x)
        {
            fmprb_t t;
            fmprb_init(t);

            fmprb_mul(t, a, b, prec);
            fmprb_mul_2exp_si(t, t, 1);
            fmprb_mul(e, a, a, prec);
            fmprb_mul(f, b, b, prec);
            fmprb_sub(e, e, f, prec);
            fmprb_swap(f, t);

            fmprb_clear(t);
        }
        else
        {
            fmprb_mul(e, a, a, prec);
            fmprb_mul(f, b, b, prec);
            fmprb_sub(e, e, f, prec);
            fmprb_mul(f, a, b, prec);
            fmprb_mul_2exp_si(f, f, 1);
        }
    }
    else
    {
        /* aliasing */
        if (z == x)
        {
            fmprb_t t, u;

            fmprb_init(t);
            fmprb_init(u);

            fmprb_mul(t, a, c, prec);
            fmprb_mul(u, a, d, prec);

            fmprb_mul(e, b, d, prec);
            fmprb_sub(e, t, e, prec);

            fmprb_mul(f, b, c, prec);
            fmprb_add(f, u, f, prec);

            fmprb_clear(t);
            fmprb_clear(u);
        }
        else if (z == y)
        {
            fmprb_t t, u;

            fmprb_init(t);
            fmprb_init(u);

            fmprb_mul(t, c, a, prec);
            fmprb_mul(u, c, b, prec);

            fmprb_mul(e, d, b, prec);
            fmprb_sub(e, t, e, prec);

            fmprb_mul(f, d, a, prec);
            fmprb_add(f, u, f, prec);

            fmprb_clear(t);
            fmprb_clear(u);
        }
        else
        {
            fmprb_t t;
            fmprb_init(t);

            fmprb_mul(e, a, c, prec);
            fmprb_mul(t, b, d, prec);
            fmprb_sub(e, e, t, prec);

            fmprb_mul(f, a, d, prec);
            fmprb_mul(t, b, c, prec);
            fmprb_add(f, f, t, prec);

            fmprb_clear(t);
        }
    }

#undef a
#undef b
#undef c
#undef d
#undef e
#undef f
}
