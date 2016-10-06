/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

static void
chi_vec_evenpart(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)
{

    ulong c3, c4, x;

    if (G->neven >= 1 && (c3 = chi->expo[0]))
    {
        for (x = 3; x < nv; x += 4)
            v[x] = c3;
    }
    if (G->neven == 2 && (c4 = chi->expo[1]))
    {
        ulong g, vx, xp;
        nmod_t pe;
        vx = c4;
        pe = G->P[1].pe;
        g = G->P[1].g;

        for (x = g; x > 1;)
        {

            for (xp = x; xp < nv; xp += pe.n)
                v[xp] = nmod_add(v[xp], vx, chi->order);

            for (xp = pe.n - x; xp < nv; xp += pe.n)
                v[xp] = nmod_add(v[xp], vx, chi->order);

            x = nmod_mul(x, g, pe);
            vx = nmod_add(vx,  c4, chi->order);
        }
    }
}

/* loop over primary components */
void
dirichlet_ui_chi_vec_primeloop(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)
{
    slong k, l;

    for (k = 0; k < nv; k++)
        v[k] = 0;

    if (G->neven)
        chi_vec_evenpart(v, G, chi, nv);

    for (l = G->neven; l < G->num; l++)
    {
        dirichlet_prime_group_struct P = G->P[l];

        /* FIXME: there may be some precomputed dlog in P if needed */
        if (P.dlog == NULL)
            dlog_vec_add(v, nv, P.g, chi->expo[l], P.pe, P.phi, chi->order);
        else
            dlog_vec_add_precomp(v, nv, P.dlog, P.g, chi->expo[l], P.pe, P.phi, chi->order);

    }
    dirichlet_ui_vec_set_null(v, G, nv);
}
