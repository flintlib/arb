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

#include "arb.h"

void
arb_agm(arb_t z, const arb_t x, const arb_t y, long prec)
{
    arb_t t, u, v, w;

    if (arb_contains_negative(x) || arb_contains_negative(y))
    {
        arb_indeterminate(z);
        return;
    }

    if (arb_is_zero(x) || arb_is_zero(y))
    {
        arb_zero(z);
        return;
    }

    arb_init(t);
    arb_init(u);
    arb_init(v);
    arb_init(w);

    arb_set(t, x);
    arb_set(u, y);

    while (!arb_overlaps(t, u) &&
            !arb_contains_nonpositive(t) &&
            !arb_contains_nonpositive(u))
    {
        arb_add(v, t, u, prec);
        arb_mul_2exp_si(v, v, -1);

        arb_mul(w, t, u, prec);
        arb_sqrt(w, w, prec);

        arb_swap(v, t);
        arb_swap(w, u);
    }

    if (!arb_is_finite(t) || !arb_is_finite(u))
    {
        arb_indeterminate(z);
    }
    else
    {
        arb_union(z, t, u, prec);
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    arb_clear(w);
}

