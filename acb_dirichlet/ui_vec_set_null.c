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
acb_dirichlet_ui_vec_set_null(ulong *v, const acb_dirichlet_group_t G, slong nv)
{
    slong k, l;
    if (G->q_even > 1)
    {
        for (k = 2; k < nv; k += 2)
            v[k] = -1;
    }

    for (l = G->neven; l < G->num; l++)
    {
        ulong p = G->P[l].p;

        for (k = p; k < nv; k += p)
            v[k] = ACB_DIRICHLET_CHI_NULL;
    }
}
