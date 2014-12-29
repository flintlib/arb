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
#include "acb.h"

void
arb_hurwitz_zeta(arb_t res, const arb_t s, const arb_t z, long prec)
{
    if (!arb_contains_si(s, 1) &&
        (arb_is_positive(z) ||
            (arb_is_int(z) && arb_is_int(s) && arb_is_nonpositive(s))))
    {
        acb_t a, b, c;

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_set_arb(a, s);
        acb_set_arb(b, z);
        acb_hurwitz_zeta(c, a, b, prec);
        arb_set(res, acb_realref(c));

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }
    else
    {
        arb_indeterminate(res);
    }
}

