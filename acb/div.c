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
    if (arb_is_zero(acb_imagref(y)))
    {
        if (arb_is_zero(acb_imagref(x)))
        {
            arb_div(acb_realref(z), acb_realref(x), acb_realref(y), prec);
            arb_zero(acb_imagref(z));
        }
        else if (arb_is_zero(acb_realref(x)))
        {
            arb_div(acb_imagref(z), acb_imagref(x), acb_realref(y), prec);
            arb_zero(acb_realref(z));
        }
        else if (z != y)
        {
            arb_div(acb_realref(z), acb_realref(x), acb_realref(y), prec);
            arb_div(acb_imagref(z), acb_imagref(x), acb_realref(y), prec);
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_set(t, acb_realref(y));
            arb_div(acb_realref(z), acb_realref(x), t, prec);
            arb_div(acb_imagref(z), acb_imagref(x), t, prec);
            arb_clear(t);
        }
    }
    else if (arb_is_zero(acb_realref(y)))
    {
        if (arb_is_zero(acb_imagref(x)))
        {
            arb_div(acb_imagref(z), acb_realref(x), acb_imagref(y), prec);
            arb_neg(acb_imagref(z), acb_imagref(z));
            arb_zero(acb_realref(z));
        }
        else if (arb_is_zero(acb_realref(x)))
        {
            arb_div(acb_realref(z), acb_imagref(x), acb_imagref(y), prec);
            arb_zero(acb_imagref(z));
        }
        else if (z != y)
        {
            arb_div(acb_realref(z), acb_realref(x), acb_imagref(y), prec);
            arb_div(acb_imagref(z), acb_imagref(x), acb_imagref(y), prec);
            arb_swap(acb_realref(z), acb_imagref(z));
            arb_neg(acb_imagref(z), acb_imagref(z));
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_set(t, acb_imagref(y));
            arb_div(acb_realref(z), acb_realref(x), t, prec);
            arb_div(acb_imagref(z), acb_imagref(x), t, prec);
            arb_swap(acb_realref(z), acb_imagref(z));
            arb_neg(acb_imagref(z), acb_imagref(z));
            arb_clear(t);
        }
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

