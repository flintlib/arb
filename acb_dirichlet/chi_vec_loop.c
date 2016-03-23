/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "acb_dirichlet.h"

/* loop over whole group */
void
acb_dirichlet_chi_vec_loop(ulong *v, ulong nv, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi)
{
    int j;
    ulong t, k;
    acb_dirichlet_conrey_t x;
    acb_dirichlet_conrey_init(x, G);
    acb_dirichlet_conrey_one(x, G);

    for (k = 0; k < nv; k++)
        v[k] = CHI_NULL;

    t = v[1] = 0;
    while ( (j = acb_dirichlet_conrey_next(x, G)) < G->num )
    {
        /* exponents were modified up to j */
        for (k = 0; k <= j; k++)
            t = (t + chi->expo[k]) % chi->order;
        if (x->n < nv)
            v[x->n] = t;
    }
    /* fix result outside primes */
    /*acb_dirichlet_vec_set_null(v, nv, G);*/
    /* copy outside modulus */
    for (k = G->q; k < nv ; k++ )
        v[k] = v[k - G->q];
    acb_dirichlet_conrey_clear(x);
}
