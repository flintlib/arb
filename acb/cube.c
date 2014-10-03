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
acb_cube(acb_t r, const acb_t z, long prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)
    if (arb_is_zero(b))
    {
        arb_pow_ui(acb_realref(r), a, 3, prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(a))
    {
        arb_pow_ui(acb_imagref(r), b, 3, prec);
        arb_neg(acb_imagref(r), acb_imagref(r));
        arb_zero(acb_realref(r));
    }
    else
    {
        arb_t t, u, v;

        arb_init(t);
        arb_init(u);
        arb_init(v);

        arb_mul(t, a, a, prec);
        arb_mul(u, b, b, prec);
        arb_set(v, t);

        /* t = a^2 - 3b^2 */
        arb_submul_ui(t, u, 3, prec);

        /* u = -(b^2 - 3a^2) */
        arb_submul_ui(u, v, 3, prec);
        arb_neg(u, u);

        arb_mul(acb_realref(r), t, a, prec);
        arb_mul(acb_imagref(r), u, b, prec);

        arb_clear(t);
        arb_clear(u);
        arb_clear(v);
    }
#undef a
#undef b
}

