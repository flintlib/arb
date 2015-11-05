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

/* [(a + b) +/- (b - a)] / 2 */

void
arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, slong prec)
{
    arf_t t;
    int inexact;

    if (arf_is_inf(a) && arf_equal(a, b))
    {
        /* [-inf, -inf] or [+inf, +inf] */
        arf_set(arb_midref(x), a);
        mag_zero(arb_radref(x));
        return;
    }

    arf_init(t);
    arf_sub(t, b, a, MAG_BITS, ARF_RND_UP);

    if (arf_sgn(t) < 0)
    {
        printf("exception: arb_set_interval_arf: endpoints not ordered\n");
        abort();
    }

    arf_get_mag(arb_radref(x), t);

    inexact = arf_add(arb_midref(x), a, b, prec, ARB_RND);
    if (inexact)
        arf_mag_add_ulp(arb_radref(x), arb_radref(x), arb_midref(x), prec);

    arb_mul_2exp_si(x, x, -1);

    arf_clear(t);
}

