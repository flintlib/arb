/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

int
acb_dirichlet_conrey_next(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G)
{
    /* update index */
    int k;
    for (k = G->num - 1; k >= 0; k--)
    {
        x->n = nmod_mul(x->n, G->generators[k], G->mod);
        x->log[k] += 1;
        if (x->log[k] < G->P[k].phi)
            break;
        x->log[k] = 0;
    }
    /* return last index modified */
    return k;
}
