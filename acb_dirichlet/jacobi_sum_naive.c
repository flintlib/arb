/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_jacobi_sum_naive(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec)
{

    ulong k1, k2, m1, m2, g, e;
    ulong * v1, * v2;
    slong *v;
    nmod_t order;
    acb_t z;

    v1 = flint_malloc(G->q * sizeof(ulong));
    v2 = flint_malloc(G->q * sizeof(ulong));

    acb_dirichlet_ui_vec_set_null(v1, G, G->q);
    acb_dirichlet_ui_chi_vec_loop(v1, G, chi1, G->q);

    acb_dirichlet_ui_vec_set_null(v2, G, G->q);
    acb_dirichlet_ui_chi_vec_loop(v2, G, chi2, G->q);

    m1 = chi1->order.n;
    m2 = chi2->order.n;
    g = n_gcd(m1, m2);
    nmod_init(&order, m1 * (m2 / g));
    m1 = order.n / m1;
    m2 = order.n / m2;

    v = flint_malloc(order.n * sizeof(slong));

    for (e = 0; e < order.n; e++)
        v[e] = 0;

    for (k1 = 2, k2 = G->q - 1; k2 > 1; k1++, k2--)
    {
        if (v1[k1] == ACB_DIRICHLET_CHI_NULL ||
            v2[k2] == ACB_DIRICHLET_CHI_NULL)
            continue;
        e = nmod_add(v1[k1] * m1, v2[k2] * m2, order);
        v[e]++;
    }

    acb_init(z);
    acb_dirichlet_nth_root(z, order.n, prec);
    acb_dirichlet_si_poly_evaluate(res, v, order.n, z, prec);

    acb_clear(z);
    flint_free(v);
    flint_free(v2);
    flint_free(v1);
}
