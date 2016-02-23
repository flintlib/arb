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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "acb_dirichlet.h"

void
acb_dirichlet_l(acb_t res, const acb_t s,
    const acb_dirichlet_group_t G, ulong m, slong prec)
{
    acb_t chi, t, u, a;
    ulong k;

    acb_init(chi);
    acb_init(t);
    acb_init(u);
    acb_init(a);

    acb_zero(t);

    for (k = 1; k <= G->q; k++)
    {
        acb_dirichlet_chi(chi, G, m, k, prec);

        if (!acb_is_zero(chi))
        {
            acb_set_ui(a, k);
            acb_div_ui(a, a, G->q, prec);
            acb_hurwitz_zeta(u, s, a, prec);
            acb_addmul(t, chi, u, prec);
        }
    }

    acb_set_ui(u, G->q);
    acb_neg(a, s);
    acb_pow(u, u, a, prec);
    acb_mul(res, t, u, prec);

    acb_clear(chi);
    acb_clear(t);
    acb_clear(u);
    acb_clear(a);
}

