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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

void
_acb_poly_taylor_shift_divconquer(acb_ptr poly, const acb_t c, slong len, slong prec)
{
    acb_struct d[2];

    if (len <= 1 || acb_is_zero(c))
        return;

    if (len == 2)
    {
        acb_addmul(poly, poly + 1, c, prec);
        return;
    }

    d[0] = *c;

    acb_init(d + 1);
    acb_one(d + 1); /* no need to free */

    _acb_poly_compose_divconquer(poly, poly, len, d, 2, prec);
}

void
acb_poly_taylor_shift_divconquer(acb_poly_t g, const acb_poly_t f,
    const acb_t c, slong prec)
{
    if (f != g)
        acb_poly_set_round(g, f, prec);

    _acb_poly_taylor_shift_divconquer(g->coeffs, c, g->length, prec);
}

