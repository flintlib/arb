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

static void
gauss_sum_non_primitive(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{

    slong k, mu = 1;
    ulong NN0 = G->q / chi->conductor;

    /* G(chi) = mu(N/N0)chi0(N/N0)G(chi0) */

    if (NN0 % 4 == 0)
    {
        acb_zero(res);
        return;
    }

    for (k = 0; k < G->num; k++)
    {
        ulong p = G->primes[k];

        if (G->exponents[k] > 1 && NN0 % (p*p) == 0)
        {
            acb_zero(res);
            return;
        }

        if (NN0 % p == 0)
            mu *= -1;
    }

    if (chi->x->n == 1)
    {
        acb_set_si(res, mu);
    }
    else
    {

        acb_dirichlet_group_t G0;
        acb_dirichlet_char_t chi0;
        acb_t z;

        /* TODO: implement efficient subgroup */
        acb_dirichlet_group_init(G0, chi->conductor);
        acb_dirichlet_char_init(chi0, G);
        acb_dirichlet_char_primitive(chi0, G0, G, chi);

        acb_init(z);
        acb_dirichlet_gauss_sum(z, G0, chi0, prec);

        acb_dirichlet_chi(res, G0, chi0, NN0, prec);

        acb_mul(res, res, z, prec);
        acb_mul_si(res, res, mu, prec);

        acb_dirichlet_group_clear(G0);
        acb_dirichlet_char_clear(chi0);
        acb_clear(z);
    }
}

void
acb_dirichlet_gauss_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{
    if (chi->conductor != G->q)
    {
        gauss_sum_non_primitive(res, G, chi, prec);
    }
    else if (chi->order <= 2)
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
