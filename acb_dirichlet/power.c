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

/* assume won't be modified */
void
acb_dirichlet_power(acb_t z, const acb_dirichlet_powers_t t, ulong n, slong prec)
{
    if (n < t->m)
    {
        /* acb_set(z, t->z + n); */
        *z = *(t->z + n);
    }
    else
    {
        ulong q, r;
        q = n / t->m;
        r = n % t->m;
        if (q >= t->M)
        {
            flint_printf("acb_dirichlet_power: power %wu not available "
                    "in table of size %wu * %wu.", n, t->m, t->M);
            abort();
        }
        acb_mul(z, t->Z + q, t->z + r, prec);
    }
}
