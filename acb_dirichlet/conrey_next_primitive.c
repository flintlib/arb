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

int
acb_dirichlet_conrey_next_primitive(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G)
{
    /* update index avoiding multiples of p except for first component
       if 8|q */
    ulong k = 0;
    if (G->neven == 2)
    {
        x->n = nmod_mul(x->n, G->generators[0], G->mod);
        x->log[0]++;
        if (x->log[0] == 1)
            return 0;
        x->log[0] = 0;
        k = 1;
    }
    for (; k < G->num ; k++)
    {
        x->n = nmod_mul(x->n, G->generators[k], G->mod);
        x->log[k]++;
        if (x->log[k] % G->P[k].p == 0)
        {
            x->n = nmod_mul(x->n, G->generators[k], G->mod);
            x->log[k]++;
        }
        if (x->log[k] < G->P[k].phi)
            break;
        x->log[k] =  1;
    }
    /* return last index modified */
    return k;
}
