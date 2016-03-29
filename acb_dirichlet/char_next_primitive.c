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
acb_dirichlet_char_next_primitive(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
    ulong k;

    acb_dirichlet_char_denormalize(chi, G);

    /* update index */
    k = 0;
    if (G->neven == 2)
    {
        chi->n = nmod_mul(chi->n, G->generators[0], G->mod);
        chi->expo[0]++;
        if (chi->expo[0] < G->expo)
            return 0;
        chi->expo[0] = 0;
        k = 1;
    }

    for (; k < G->num ; k++)
    {
        chi->n = nmod_mul(chi->n, G->generators[k], G->mod);
        chi->expo[k] += G->PHI[k];
        if (chi->expo[k] % G->primes[k] == 0)
        {
            chi->n = nmod_mul(chi->n, G->generators[k], G->mod);
            chi->expo[k] += G->PHI[k];
        }
        if (chi->expo[k] < G->expo)
            break;
        chi->expo[k] = G->PHI[k];  
    }

    acb_dirichlet_char_normalize(chi, G);

    /* return last index modified */
    return k;
}
