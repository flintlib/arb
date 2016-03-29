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
acb_dirichlet_pairing(const acb_dirichlet_group_t G, ulong m, ulong n)
{
    ulong x;
    acb_dirichlet_conrey_t a, b;

    if (n_gcd(G->q, m) > 1 || n_gcd(G->q, n) > 1)
        return ACB_DIRICHLET_CHI_NULL;

    acb_dirichlet_conrey_init(a, G);
    acb_dirichlet_conrey_init(b, G);
    acb_dirichlet_conrey_log(a, G, m);
    acb_dirichlet_conrey_log(b, G, n);

    x = acb_dirichlet_pairing_conrey(G, a, b);

    acb_dirichlet_conrey_clear(a);
    acb_dirichlet_conrey_clear(b);

    return x;
}
