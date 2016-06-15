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

/* order of an element knowing the factorization of a multiple */
ulong
nmod_order_precomp(ulong a, nmod_t mod, ulong expo, n_factor_t fac)
{
    int k;
    ulong pe, ap, order = 1;
    for (k = 0; k < fac.num; k++)
    {
        pe = n_pow(fac.p[k], fac.exp[k]);
        ap = nmod_pow_ui(a, expo / pe, mod);
        while ( ap != 1)
        {
            ap = nmod_pow_ui(ap, fac.p[k], mod);
            order *= fac.p[k];
        }
    }
    return order;
}

ulong
acb_dirichlet_ui_order(const acb_dirichlet_group_t G, ulong a)
{
    n_factor_t fac;

    n_factor_init(&fac);
    n_factor(&fac, G->expo, 1);
    return nmod_order_precomp(a, G->mod, G->expo, fac);
}
