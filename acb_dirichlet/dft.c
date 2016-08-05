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

/* dft, lexicographic conrey indexing, array size G->phi_q */
void
acb_dirichlet_dft_conrey(acb_ptr w, acb_srcptr v, const acb_dirichlet_group_t G, slong prec)
{
    slong k, l, * cyc;
    cyc = flint_malloc(G->num * sizeof(slong));
    for (k = 0, l = G->num - 1; l >= 0; k++, l--)
        cyc[k] = G->P[k].phi;

    acb_dirichlet_dft_prod(w, v, cyc, G->num, prec);
    flint_free(cyc);
}

/* dft, number indexing, array size G->q */
void
acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const acb_dirichlet_group_t G, slong prec)
{
    ulong i, len;
    acb_ptr t1, t2;
    acb_dirichlet_conrey_t x;

    len = G->phi_q;
    t1 = flint_malloc(len * sizeof(acb_struct));

    acb_dirichlet_conrey_init(x, G);
    acb_dirichlet_conrey_one(x, G);
    for (i = 0; i < len; i++)
    {
        t1[i] = v[x->n];
        acb_dirichlet_conrey_next(x, G);
    };

    t2 = _acb_vec_init(len);
    acb_dirichlet_dft_conrey(t2, t1, G, prec);

    acb_dirichlet_conrey_one(x, G);
    for (i = 0; i < len; i++)
    {
        acb_set(w + x->n, t2 + i);
        acb_dirichlet_conrey_next(x, G);
    };

    _acb_vec_clear(t2, len);
    acb_dirichlet_conrey_clear(x);
    flint_free(t1);
}
