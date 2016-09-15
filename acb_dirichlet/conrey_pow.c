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
acb_dirichlet_conrey_pow(acb_dirichlet_conrey_t c, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t a, ulong n)
{
    ulong k;
    for (k = 0; k < G->num ; k++)
        c->log[k] = n_mulmod2(a->log[k], n, G->P[k].phi);
    c->n = nmod_pow_ui(a->n, n, G->mod);
}
