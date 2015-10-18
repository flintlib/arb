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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

void
acb_hypgeom_2f1_corner(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, long prec)
{
    acb_t aa, bb, cc, z1, z2, f1, f2;
    int upper;

    acb_init(aa); acb_init(bb); acb_init(cc);
    acb_init(z1); acb_init(z2); acb_init(f1); acb_init(f2);

    acb_add_ui(aa, a, 1, prec);
    acb_add_ui(bb, b, 1, prec);
    acb_add_ui(cc, c, 1, prec);

    upper = arb_is_positive(acb_imagref(z));

    /* 0 -> 0.5 +/- 0.5i -> 0.5 +/- 0.75i -> z */
    acb_set_d_d(z1, 0.5, upper ? 0.5 : -0.5);
    acb_set_d_d(z2, 0.5, upper ? 0.75 : -0.75);

    acb_hypgeom_2f1(f1, a, b, c, z1, regularized, prec);
    acb_hypgeom_2f1(f2, aa, bb, cc, z1, regularized, prec);
    acb_mul(f2, f2, a, prec);
    acb_mul(f2, f2, b, prec);
    if (!regularized)
        acb_div(f2, f2, c, prec);

    acb_hypgeom_2f1_continuation(f1, f2, a, b, c, z1, z2, f1, f2, prec);

    acb_set(z1, z2);
    acb_set(z2, z);

    acb_hypgeom_2f1_continuation(f1, f2, a, b, c, z1, z2, f1, f2, prec);

    acb_set(res, f1);

    acb_clear(aa); acb_clear(bb); acb_clear(cc);
    acb_clear(z1); acb_clear(z2); acb_clear(f1); acb_clear(f2);
}

