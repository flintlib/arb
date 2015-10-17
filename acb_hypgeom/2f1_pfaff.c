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
acb_hypgeom_2f1_pfaff(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, long prec)
{
    acb_t t, u, v;

    acb_init(t);
    acb_init(u);
    acb_init(v);

    acb_sub_ui(t, z, 1, prec); /* t = z-1 */
    acb_div(u, z, t, prec); /* u = z/(z-1) */
    acb_neg(t, t);
    acb_neg(v, a);
    acb_pow(t, t, v, prec); /* t = (1-z)^-a */
    acb_sub(v, c, b, prec); /* v = c-b */

    acb_hypgeom_2f1_direct(res, a, v, c, u, regularized, prec);
    acb_mul(res, res, t, prec);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

