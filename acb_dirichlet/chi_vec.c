/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_chi_vec(acb_ptr v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv, slong prec)
{
    slong k;
    ulong * a;
    acb_dirichlet_powers_t t;

    a = flint_malloc(nv * sizeof(ulong));
    acb_dirichlet_ui_chi_vec(a, G, chi, nv);

    acb_dirichlet_powers_init(t, chi->order.n, nv, prec);

    acb_zero(v + 0);
    for (k = 0; k < nv; k++)
    {
        if (a[k] != ACB_DIRICHLET_CHI_NULL)
            acb_dirichlet_power(v + k, t, a[k], prec);
        else
            *(v + k) = *(v + 0);
    }

    acb_dirichlet_powers_clear(t);
    flint_free(a);
}
