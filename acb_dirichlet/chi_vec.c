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
acb_dirichlet_chi_vec(acb_ptr v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv, slong prec)
{
    slong k;
    ulong * a;
    acb_dirichlet_powers_t t;

    a = flint_malloc(nv * sizeof(ulong));
    acb_dirichlet_ui_chi_vec(a, G, chi, nv);

    acb_dirichlet_powers_init(t, chi->order.n, nv, prec);

    acb_zero(v + 0);
    for (k = 0; k < nv; k++)
    {
        if (a[k] != ACB_DIRICHLET_CHI_NULL)
            acb_dirichlet_power(v + k, t, a[k], prec);
        else
            *(v + k) = *(v + 0);
    }

    acb_dirichlet_powers_clear(t);
    flint_free(a);
}
