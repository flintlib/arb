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
acb_dirichlet_char_set_expo(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
    slong k;
    for (k = 0; k < G->num; k++)
        chi->expo[k] = (chi->x->log[k] * G->PHI[k]) % G->expo;
}

void
acb_dirichlet_char_normalize(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
    ulong k, g;
    g = G->expo;

    for (k = 0; k < G->num; k++)
        g = n_gcd(g, chi->expo[k]);

    for (k = 0; k < G->num; k++)
        chi->expo[k] = chi->expo[k] / g;

    chi->order = G->expo / g;
}

void
acb_dirichlet_char_denormalize(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
    ulong k, g;
    g = G->expo / chi->order;

    for (k = 0; k < G->num; k++)
        chi->expo[k] *= g;
}
