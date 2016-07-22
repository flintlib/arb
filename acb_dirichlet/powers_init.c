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
_acb_dirichlet_powers_init(acb_dirichlet_powers_t t, ulong order, slong size, slong depth, slong prec)
{
    slong k;
    t->order = order;
    acb_init(t->z);
    acb_dirichlet_nth_root(t->z, order, prec);
    t->size = size;
    t->depth = depth;
    if (depth)
    {
        acb_struct * z;
        z = t->z + 0;
        t->Z = flint_malloc(depth * sizeof(acb_ptr));
        for (k = 0; k < depth; k++)
        {
            t->Z[k] = _acb_vec_init(size);
            _acb_vec_set_powers(t->Z[k], z, size, prec);
            z = t->Z[k] + 1;
        }
    }
    else
    {
        t->Z = NULL;
    }
}

void
acb_dirichlet_powers_init(acb_dirichlet_powers_t t, ulong order, slong num, slong prec)
{
    slong size, depth;

    depth = (num > order) ? 1 : n_flog(order, num);
    size = n_root(order, depth) + 1;

    _acb_dirichlet_powers_init(t, order, size, depth, prec);
}
