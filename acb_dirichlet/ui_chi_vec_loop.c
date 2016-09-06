/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* loop over whole group */
void
acb_dirichlet_ui_chi_vec_loop(ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv)
{
    int j;
    ulong t;
    slong k;
    acb_dirichlet_conrey_t x;
    acb_dirichlet_conrey_init(x, G);
    acb_dirichlet_conrey_one(x, G);

    for (k = 0; k < nv; k++)
        v[k] = ACB_DIRICHLET_CHI_NULL;

    t = v[1] = 0;

    while ( (j = acb_dirichlet_conrey_next(x, G)) >= 0 )
    {
        /* exponents were modified up to j */
        for (k = G->num - 1; k >= j; k--)
            t = nmod_add(t, chi->expo[k], chi->order);

        if (x->n < nv)
            v[x->n] = t;
    }

    /* fix result outside primes */
    /*acb_dirichlet_vec_set_null(v, G, nv);*/
    /* copy outside modulus */

    for (k = G->q; k < nv ; k++ )
        v[k] = v[k - G->q];

    acb_dirichlet_conrey_clear(x);
}
