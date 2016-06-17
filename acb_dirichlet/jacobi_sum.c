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
acb_dirichlet_jacobi_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec)
{
    if (chi1->x->n == 1 || chi2->x->n == 1)
    {
        if (chi1->x->n == 1 && chi2->x->n == 1)
            acb_set_si(res, G->q - 2);
        else
            acb_set_si(res, - 1);
    }
    else if (nmod_mul(chi1->x->n, chi2->x->n, G->mod) == 1)
    {
        if (chi1->parity)
            acb_one(res);
        else
            acb_set_si(res, -1);
    }
    else
    {
        acb_dirichlet_char_t chi12;
        acb_t tmp;

        acb_dirichlet_char_init(chi12, G);
        acb_init(tmp);

        acb_dirichlet_gauss_sum(res, G, chi1, prec);
        acb_dirichlet_gauss_sum(tmp, G, chi2, prec);
        acb_mul(res, res, tmp, prec);
        acb_dirichlet_char_mul(chi12, G, chi1, chi2);
        acb_dirichlet_gauss_sum(tmp, G, chi12, prec);
        acb_div(res, res, tmp, prec);

        acb_dirichlet_char_clear(chi12);
        acb_clear(tmp);
    }

}
