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
dirichlet_char_set_expo(dirichlet_char_t chi, const dirichlet_group_t G)
{
    slong k;
    for (k = 0; k < G->num; k++)
        /* no overflow: log[k] < phi[k] and G->expo = phi[k] * PHI[k] */
        chi->expo[k] = chi->x->log[k] * G->PHI[k];
}

void
dirichlet_char_normalize(dirichlet_char_t chi, const dirichlet_group_t G)
{
    ulong k, g;
    g = G->expo;

    for (k = 0; k < G->num; k++)
        g = n_gcd(g, chi->expo[k]);

    for (k = 0; k < G->num; k++)
        chi->expo[k] = chi->expo[k] / g;

    nmod_init(&chi->order, G->expo / g);
}

void
dirichlet_char_denormalize(dirichlet_char_t chi, const dirichlet_group_t G)
{
    ulong k, g;
    g = G->expo / chi->order.n;

    for (k = 0; k < G->num; k++)
        chi->expo[k] *= g;
}
