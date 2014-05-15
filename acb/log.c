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

#include "acb.h"

void
acb_log(acb_t r, const acb_t z, long prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)

    arb_t t, u;

    arb_init(t);
    arb_init(u);

    arb_mul(t, a, a, prec);
    arb_mul(u, b, b, prec);
    arb_add(t, t, u, prec);

    if (arb_contains_zero(t) || arf_sgn(arb_midref(t)) < 0)
    {
        arb_zero_pm_inf(t);
    }
    else
    {
        arb_log(t, t, prec);
    }

    acb_arg(u, z, prec);

    arb_mul_2exp_si(acb_realref(r), t, -1);
    arb_set(acb_imagref(r), u);

    arb_clear(t);
    arb_clear(u);
}

