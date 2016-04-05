/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_l(acb_t res, const acb_t s,
    const acb_dirichlet_group_t G, ulong m, slong prec)
{
    acb_t chi, t, u, a;
    acb_dirichlet_conrey_t cm, cn;

    acb_init(chi);
    acb_dirichlet_conrey_init(cm, G);
    acb_dirichlet_conrey_init(cn, G);
    acb_init(t);
    acb_init(u);
    acb_init(a);

    acb_dirichlet_conrey_log(cm, G, m);
    acb_dirichlet_conrey_one(cn, G);
    acb_zero(t);

    while (1) {
        /* todo: use n_dirichlet_chi and precomputed roots instead */
        acb_dirichlet_pairing_conrey(chi, G, cm, cn, prec);

        acb_set_ui(a, cn->n);
        acb_div_ui(a, a, G->q, prec);
        acb_hurwitz_zeta(u, s, a, prec);
        acb_addmul(t, chi, u, prec);

        if (acb_dirichlet_conrey_next(cn, G) == G->num)
          break;
    }

    acb_set_ui(u, G->q);
    acb_neg(a, s);
    acb_pow(u, u, a, prec);
    acb_mul(res, t, u, prec);

    acb_dirichlet_conrey_clear(cm);
    acb_dirichlet_conrey_clear(cn);

    acb_clear(chi);
    acb_clear(t);
    acb_clear(u);
    acb_clear(a);
}
