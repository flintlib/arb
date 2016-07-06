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

ulong
acb_dirichlet_conrey_conductor(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)
{
    int k, f;
    ulong cond = 1;

    if (G->neven >= 1 && x->log[0] == 1)
        cond = 4;
    if (G->neven == 2 && x->log[1])
    {
        ulong l2 = x->log[1];
        f = n_remove(&l2, 2);
        cond = 1 << (G->P[1].e - f);
    }

    for (k = G->neven; k < G->num; k++)
    {
        if (x->log[k])
        {
            ulong p, lp;
            p = G->P[k].p;
            lp = x->log[k];
            f = n_remove(&lp, p);
            if (f)
                cond *= n_pow(p, G->P[k].e - f);
            else
                cond *= G->P[k].pe.n;
        }
    }

    return cond;
}
