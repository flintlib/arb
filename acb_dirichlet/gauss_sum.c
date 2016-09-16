/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static void
gauss_sum_non_primitive(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{

    slong k, mu = 1;
    ulong NN0 = G->q / chi->conductor;

    /* G(chi) = mu(N/N0)chi0(N/N0)G(chi0) */

    if (NN0 % 2 == 0)
    {
        if (G->q % 4)
            mu = -1;
        else
        {
            acb_zero(res);
            return;
        }
    }

    for (k = G->neven; k < G->num; k++)
    {
        ulong p = G->P[k].p;

        if (G->P[k].e > 1 && NN0 % (p*p) == 0)
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

        acb_dirichlet_subgroup_init(G0, G, chi->conductor);
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
    /* TODO: no need, factor also does it... */
    if (chi->conductor != G->q)
    {
        gauss_sum_non_primitive(res, G, chi, prec);
    }
    else if (chi->order.n <= 2)
    {
        acb_dirichlet_gauss_sum_order2(res, chi, prec);
    }
    else if (G->num  > 1 && G->num > G->neven)
    {
        acb_dirichlet_gauss_sum_factor(res, G, chi, prec);
    }
    else
    {
        if (acb_dirichlet_theta_length_d(G->q, 1, prec) > G->q)
            acb_dirichlet_gauss_sum_naive(res, G, chi, prec);
        else
            acb_dirichlet_gauss_sum_theta(res, G, chi, prec);
    }
}

void
acb_dirichlet_gauss_sum_ui(acb_t res, const acb_dirichlet_group_t G, ulong a, slong prec)
{
    acb_dirichlet_char_t chi;
    acb_dirichlet_char_init(chi, G);
    acb_dirichlet_char(chi, G, a);
    acb_dirichlet_gauss_sum(res, G, chi, prec);
    acb_dirichlet_char_clear(chi);
}
