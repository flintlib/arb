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
#include "acb_poly.h"

void
acb_dirichlet_power(acb_t z, const acb_dirichlet_powers_t t, ulong n, slong prec)
{
    if (!t->depth)
    {
        acb_pow_ui(z, t->z, n, prec);
    }
    else
    {
        slong k;
        ulong r;
        r = n % t->size;
        n = n / t->size;
        acb_set(z, t->Z[0] + r);
        for (k = 1; k < t->depth && n; k++)
        {
            r = n % t->size;
            n = n / t->size;
            acb_mul(z, z, t->Z[k] + r, prec);
        }
    }
}
