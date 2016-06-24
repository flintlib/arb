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
acb_dirichlet_conrey_primitive(acb_dirichlet_conrey_t y, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x, ulong cond)
{
    int k, l, f;
    
    l = 0;
    if (cond % 4 == 0)
    {
        y->log[l++] = x->log[0];
        
        if (cond % 8 == 0)
        {
            ulong l2 = x->log[1];
            f = n_remove(&l2, 2);
            y->log[l++] = l2;
        }
    }

    for (k = G->neven; k < G->num; k++)
    {
        if (x->log[k])
        {
            ulong p, lp;
            p = G->primes[k];
            if (cond % p == 0)
            {
                lp = x->log[k];
                f = n_remove(&lp, p);
                y->log[l++] = lp;
            }
        }
    }

}
