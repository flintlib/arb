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

/* dft, use conrey indices */
void
acb_dirichlet_dft_conrey(acb_ptr w, acb_srcptr v, const acb_dirichlet_group_t G, slong prec)
{
    slong k, l, * cyc;
    cyc = flint_malloc(G->num * sizeof(slong));
    for (k = 0, l = G->num - 1; l >= 0; k++, l--)
        cyc[k] = G->P[l].phi;

    acb_dirichlet_dft_prod(w, v, cyc, G->num, prec);
    flint_free(cyc);
}
