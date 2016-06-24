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
acb_dirichlet_si_poly_evaluate(acb_t res, slong * v, slong len, const acb_t z, slong prec)
{
    slong k, r, m;
    acb_t sq;
    acb_ptr zk;

    if (len < 3)
    {
        if (len == 0)
        {
            acb_zero(res);
        }
        else if (len == 1)
        {
            acb_set_si(res, v[0]);
        }
        else if (len == 2)
        {
            acb_mul_si(res, z, v[1], prec);
            acb_add_si(res, res, v[0], prec);
        }
        return;
    }

    m = n_sqrt(len) + 1;

    zk = _acb_vec_init(m + 1);
    _acb_vec_set_powers(zk, z, m + 1, prec);

    acb_init(sq);
    acb_zero(res);

    k = len - 1;
    r = k % m;
    for (; k >= 0; r = m - 1)
    {
        acb_zero(sq);
        for (; r >= 0; r--, k--)
            acb_addmul_si(sq, zk + r, v[k], prec);
        acb_mul(res, res, zk + m, prec);
        acb_add(res, res, sq, prec);
    }

    _acb_vec_clear(zk, m + 1);
    acb_clear(sq);
}
