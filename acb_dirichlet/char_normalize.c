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
acb_dirichlet_char_set_expo(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
    slong k;
    for (k = 0; k < G->num; k++)
        chi->expo[k] = (chi->x->log[k] * G->PHI[k]) % G->expo;
}

void
acb_dirichlet_char_normalize(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
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
acb_dirichlet_char_denormalize(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
    ulong k, g;
    g = G->expo / chi->order.n;

    for (k = 0; k < G->num; k++)
        chi->expo[k] *= g;
}
