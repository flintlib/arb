/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

void
dirichlet_conrey_mul(dirichlet_conrey_t c, const dirichlet_group_t G, const dirichlet_conrey_t a, const dirichlet_conrey_t b)
{
    ulong k;
    for (k = 0; k < G->num ; k++)
        c->log[k] = (a->log[k] + b->log[k]) % G->P[k].phi;
    c->n = nmod_mul(a->n, b->n, G->mod);
}
