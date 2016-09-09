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
acb_dirichlet_conrey_primitive(acb_dirichlet_conrey_t y, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x, ulong cond)
{
    slong k, l;

    l = 0;
    if (cond % 4 == 0)
    {
        y->log[l++] = x->log[0];

        if (cond % 8 == 0)
        {
            ulong l2 = x->log[1];
            n_remove(&l2, 2);
            y->log[l++] = l2;
        }
    }

    for (k = G->neven; k < G->num; k++)
    {
        if (x->log[k])
        {
            ulong p, lp;
            p = G->P[k].p;
            if (cond % p == 0)
            {
                lp = x->log[k];
                n_remove(&lp, p);
                y->log[l++] = lp;
            }
        }
    }

}
