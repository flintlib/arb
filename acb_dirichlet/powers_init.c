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
#include "acb_poly.h"

void
acb_dirichlet_powers_init(acb_dirichlet_powers_t t, ulong order, slong num, slong prec)
{
    ulong m;
    acb_t zeta;
    t->order = order;

    m = (num == 1) ? 1 : num * (prec / 64 + 1);
    if (m > order)
        m = order;
    t->m = m;

    acb_init(zeta);
    acb_dirichlet_nth_root(zeta, order, prec);
    t->z = _acb_vec_init(m);
    _acb_vec_set_powers(t->z, zeta, m, prec);

    if (order > m)
    {
        t->M = (order / m) + 1;
        t->Z = _acb_vec_init(t->M);
        acb_pow_ui(zeta, zeta, m, prec);
        _acb_vec_set_powers(t->Z, zeta, t->M, prec);
    }
    else
    {
        t->M = 0;
        t->Z = NULL;
    }
    acb_clear(zeta);
}
