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
fmpcb_log(fmpcb_t r, const fmpcb_t z, long prec)
{
#define a fmpcb_realref(z)
#define b fmpcb_imagref(z)

    fmprb_t t, u;

    fmprb_init(t);
    fmprb_init(u);

    fmprb_mul(t, a, a, prec);
    fmprb_mul(u, b, b, prec);
    fmprb_add(t, t, u, prec);

    if (fmprb_contains_zero(t) || fmpr_sgn(fmprb_midref(t)) < 0)
    {
        fmpr_zero(fmprb_midref(t));
        fmpr_pos_inf(fmprb_radref(t));
    }
    else
    {
        fmprb_log(t, t, prec);
    }

    fmpcb_arg(u, z, prec);

    fmprb_mul_2exp_si(fmpcb_realref(r), t, -1);
    fmprb_set(fmpcb_imagref(r), u);

    fmprb_clear(t);
    fmprb_clear(u);
}

