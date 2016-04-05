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

void
chi_vec_evenpart(ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv)
{

    ulong c3, c4, x;

    if (G->neven >= 1 && (c3 = chi->expo[0]))
    {
        for (x = 3; x < nv; x += 4)
            v[x] = c3;
    }
    if (G->neven == 2 && (c4 = chi->expo[1]))
    {
        ulong g, pe, vx, xp;
        vx = c4;
        pe = G->primepowers[1];
        g = G->generators[1] % pe;

        for (x = g; x > 1;)
        {

            for (xp = x; xp < nv; xp += pe)
                v[xp] = (v[xp] + vx) % chi->order;

            for (xp = pe - x; xp < nv; xp += pe)
                v[xp] = (v[xp] + vx) % chi->order;

            x = (x * g) % pe;
            vx = (vx + c4) % chi->order;
        }
    }
}

/* loop over primary components */
void
acb_dirichlet_chi_vec_primeloop(ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv)
{
    slong k, l;

    for (k = 1; k < nv; k++)
        v[k] = 0;

    if (G->neven)
        chi_vec_evenpart(v, G, chi, nv);

    for (l = G->neven; l < G->num; l++)
    {
        ulong p, pe, g, x, vx, vp, xp;
        p = G->primes[l];
        pe = G->primepowers[l];
        g = G->generators[l] % pe;
        vx = vp = chi->expo[l];

        if (vp == 0)
            continue;
        /* for each x = g^j mod p^e,
         * set a[x] += j*vp
         * and use periodicity */
        for (x = g; x > 1;)
        {

            for (xp = x; xp < nv; xp += pe)
                v[xp] = (v[xp] + vx) % chi->order;

            x = (x * g) % pe;
            vx = (vx + vp) % chi->order;
        }
    }
    acb_dirichlet_vec_set_null(v, G, nv);
}
