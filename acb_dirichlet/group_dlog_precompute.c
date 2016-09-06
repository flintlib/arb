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
acb_dirichlet_prime_group_dlog_precompute(acb_dirichlet_prime_group_struct * P, ulong num)
{
    P->dlog = flint_malloc(sizeof(dlog_precomp_t));
    dlog_precomp_modpe_init(P->dlog, P->g, P->p, P->e, P->pe.n, num);
}


void
acb_dirichlet_group_dlog_precompute(acb_dirichlet_group_t G, ulong num)
{
    slong k;
    for (k = 0; k < G->num; k++)
        acb_dirichlet_prime_group_dlog_precompute(&G->P[k], num);
}
