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
fmpcb_inv(fmpcb_t z, const fmpcb_t x, long prec)
{
#define a fmpcb_realref(x)
#define b fmpcb_imagref(x)
#define c fmpcb_realref(z)
#define d fmpcb_imagref(z)

    if (fmprb_is_zero(b))
    {
        fmprb_inv(c, a, prec);
        fmprb_zero(d);
    }
    else if (fmprb_is_zero(a))
    {
        fmprb_inv(d, b, prec);
        fmprb_neg(d, d);
        fmprb_zero(c);
    }
    else
    {
        fmprb_t t;
        fmprb_init(t);

        fmprb_mul(t, a, a, prec);
        fmprb_addmul(t, b, b, prec);

        fmprb_div(c, a, t, prec);
        fmprb_div(d, b, t, prec);

        fmprb_neg(d, d);

        fmprb_clear(t);
    }

#undef a
#undef b
#undef c
#undef d
}

