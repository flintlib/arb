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

/* loop over primary components */
void
acb_dirichlet_chi_vec_primeloop(ulong *v, ulong nv, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi)
{
    ulong k, l;

    for (k = 1; k < nv; k++)
        v[k] = 0;

    for (l = 1; l < G->num; l++)
    {
        long p, pe, g, x, vp, xp;
        long j, vj;
        p = G->primes[l];
        pe = G->primepowers[l];
        g = G->generators[l] % pe;
        vj = vp = chi->expo[l];
        /* for each x = g^j mod p^e,
         * set a[x] += j*vp
         * and use periodicity */
        for (j = 1, x = g; x > 1; j++)
        {
            for (xp = x; xp < nv; xp += pe)
                v[xp] = (v[xp] + vj) % chi->order;

            x = (x*g) % pe;
            vj = (vj + vp) % chi->order;
        }
    }
    acb_dirichlet_vec_set_null(v, nv, G);
}
