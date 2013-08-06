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

#include "fmprb.h"

void
fmprb_agm(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    fmprb_t t, u, v, w;

    if (fmprb_contains_negative(x) || fmprb_contains_negative(y))
    {
        fmprb_indeterminate(z);
        return;
    }

    if (fmprb_is_zero(x) || fmprb_is_zero(y))
    {
        fmprb_zero(z);
        return;
    }

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);
    fmprb_init(w);

    fmprb_set(t, x);
    fmprb_set(u, y);

    while (!fmprb_overlaps(t, u) &&
            !fmprb_contains_nonpositive(t) &&
            !fmprb_contains_nonpositive(u))
    {
        fmprb_add(v, t, u, prec);
        fmprb_mul_2exp_si(v, v, -1);

        fmprb_mul(w, t, u, prec);
        fmprb_sqrt(w, w, prec);

        fmprb_swap(v, t);
        fmprb_swap(w, u);
    }

    if (!fmprb_is_finite(t) || !fmprb_is_finite(u))
    {
        fmprb_indeterminate(z);
    }
    else
    {
        fmprb_union(z, t, u, prec);
    }

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
    fmprb_clear(w);
}

