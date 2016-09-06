/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

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
