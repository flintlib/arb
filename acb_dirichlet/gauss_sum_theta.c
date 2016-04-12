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
acb_dirichlet_gauss_sum_theta(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{
    arb_t x;
    acb_t eps;

    arb_init(x);

    if ((G->q == 300 && (chi->n == 271 || chi->n == 131))
            || (G->q == 600 && (chi->n == 11 || chi->n == 91)))
    {
        /* or could use l'Hopital rule */
        acb_dirichlet_gauss_sum_naive(res, G, chi, prec);
        return;
    }

    arb_one(x);
    acb_dirichlet_chi_theta(res, G, chi, x, prec);
    acb_init(eps);
    acb_conj(eps, res);
    acb_div(eps, res, eps, prec);

    if (chi->parity)
    {
        arb_zero(acb_realref(res));
        arb_sqrt_ui(acb_imagref(res), G->q, prec);
    }
    else
    {
        arb_zero(acb_imagref(res));
        arb_sqrt_ui(acb_realref(res), G->q, prec);
    }

    acb_mul(res, res, eps, prec);

    arb_clear(x);
    acb_clear(eps);
}
