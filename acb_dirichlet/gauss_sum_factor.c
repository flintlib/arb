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
acb_dirichlet_gauss_sum_factor(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{

    slong k;
    acb_t tmp;

    for (k = (G->neven == 2); k < G->num; k++)
    {
        /* if e > 1 and not primitive, 0 */
        if (chi->x->log[k] % G->P[k].p == 0 && G->P[k].e > 1)
        {
            acb_zero(res);
            return;
        }
    }

    /* factor */
    acb_one(res);
    acb_init(tmp);

    for (k = (G->neven == 2); k < G->num; k++)
    {
        ulong pe = G->P[k].pe.n;
        acb_dirichlet_group_t Gp;
        acb_dirichlet_char_t chip;

        acb_dirichlet_subgroup_init(Gp, G, pe);
        acb_dirichlet_char_init(chip, Gp);

        chip->x->n = chi->x->n % pe;

        if (k == 1 && G->neven == 2)
        {
            chip->x->log[0] = chi->x->log[0];
            chip->x->log[1] = chi->x->log[1];
        }
        else
            chip->x->log[0] = chi->x->log[k];

        acb_dirichlet_char_conrey(chip, Gp, NULL);

        /*  chi_pe(a, q/pe) * G_pe(a) */
        acb_dirichlet_gauss_sum(tmp, Gp, chip, prec);
        acb_mul(res, res, tmp, prec);

        acb_dirichlet_chi(tmp, Gp, chip, (G->q / pe) % pe, prec);
        acb_mul(res, res, tmp, prec);

        acb_dirichlet_char_clear(chip);
        acb_dirichlet_group_clear(Gp);
    }

    if (G->q_even == 2)
        acb_neg(res, res);

    acb_clear(tmp);
}
