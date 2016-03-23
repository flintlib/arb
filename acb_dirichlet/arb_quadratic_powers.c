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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "acb_dirichlet.h"

void
acb_dirichlet_arb_quadratic_powers(arb_ptr v, slong nv, const arb_t x, slong prec)
{
    slong i;
    arb_t dx, x2;
    arb_init(dx);
    arb_init(x2);
    arb_set(dx, x);
    arb_mul(x2, x, x, prec);
    for (i = 0; i < nv; i++)
    {
        if (i == 0)
            arb_one(v + i);
        else if (i == 1)
            arb_set_round(v + i, x, prec);
        else
        {
            arb_mul(dx, dx, x2, prec);
            arb_mul(v + i, v + i - 1, dx, prec);
        }
    }
    arb_clear(dx);
    arb_clear(x2);
}
