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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "acb.h"

void
acb_div(acb_t z, const acb_t x, const acb_t y, long prec)
{
#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(y)
#define d acb_imagref(y)
    if (arb_is_zero(d))
    {
        if (arb_is_zero(b))
        {
            arb_div(acb_realref(z), a, c, prec);
            arb_zero(acb_imagref(z));
        }
        else if (arb_is_zero(a))
        {
            arb_div(acb_imagref(z), b, c, prec);
            arb_zero(acb_realref(z));
        }
        else if (z != y)
        {
            arb_div(acb_realref(z), a, c, prec);
            arb_div(acb_imagref(z), b, c, prec);
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_set(t, c);
            arb_div(acb_realref(z), a, t, prec);
            arb_div(acb_imagref(z), b, t, prec);
            arb_clear(t);
        }
    }
    else if (arb_is_zero(c))
    {
        if (arb_is_zero(b))
        {
            arb_div(acb_imagref(z), a, d, prec);
            arb_neg(acb_imagref(z), acb_imagref(z));
            arb_zero(acb_realref(z));
        }
        else if (arb_is_zero(a))
        {
            arb_div(acb_realref(z), b, d, prec);
            arb_zero(acb_imagref(z));
        }
        else if (z != y)
        {
            arb_div(acb_realref(z), a, d, prec);
            arb_div(acb_imagref(z), b, d, prec);
            arb_swap(acb_realref(z), acb_imagref(z));
            arb_neg(acb_imagref(z), acb_imagref(z));
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_set(t, d);
            arb_div(acb_realref(z), a, t, prec);
            arb_div(acb_imagref(z), b, t, prec);
            arb_swap(acb_realref(z), acb_imagref(z));
            arb_neg(acb_imagref(z), acb_imagref(z));
            arb_clear(t);
        }
    }
    else
    {
        if (prec > 256 && acb_bits(y) <= prec / 2)
        {
            arb_t t, u, v;

            arb_init(t);
            arb_init(u);
            arb_init(v);

            arb_mul(t, c, c, prec);
            arb_addmul(t, d, d, prec);

            arb_mul(u, a, c, prec);
            arb_addmul(u, b, d, prec);

            arb_mul(v, b, c, prec);
            arb_submul(v, a, d, prec);

            arb_div(acb_realref(z), u, t, prec);
            arb_div(acb_imagref(z), v, t, prec);

            arb_clear(t);
            arb_clear(u);
            arb_clear(v);
        }
        else
        {
            acb_t t;
            acb_init(t);
            acb_inv(t, y, prec);
            acb_mul(z, x, t, prec);
            acb_clear(t);
        }
    }
#undef a
#undef b
#undef c
#undef d
}

