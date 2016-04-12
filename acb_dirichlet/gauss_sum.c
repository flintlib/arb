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
acb_dirichlet_gauss_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{
    if (chi->order <= 2)
    {
        if (chi->parity)
        {
            arb_sqrt_ui(acb_imagref(res), G->q, prec);
            arb_zero(acb_realref(res));
        }
        else
        {
            arb_sqrt_ui(acb_realref(res), G->q, prec);
            arb_zero(acb_imagref(res));
        }
    }
    else
    {
        if (acb_dirichlet_theta_length_d(G->q, 1, prec) > G->q)
            acb_dirichlet_gauss_sum_naive(res, G, chi, prec);
        else
            acb_dirichlet_gauss_sum_theta(res, G, chi, prec);
    }
}
