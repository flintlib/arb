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
acb_dirichlet_l_vec_hurwitz(acb_ptr res, const acb_t s,
    const acb_dirichlet_group_t G, slong prec)
{
    acb_t a, qs;
    acb_ptr zeta, z;
    acb_dirichlet_conrey_t cn;

    acb_dirichlet_conrey_init(cn, G);
    acb_init(qs);
    acb_init(a);

    prec += n_clog(G->phi_q, 2);

    acb_set_ui(qs, G->q);
    acb_neg(a, s);
    acb_pow(qs, qs, a, prec);

    zeta = z = _acb_vec_init(G->phi_q);
    acb_dirichlet_conrey_one(cn, G);
    do {

        acb_set_ui(a, cn->n);
        acb_div_ui(a, a, G->q, prec);
        acb_hurwitz_zeta(z, s, a, prec);
        acb_mul(z, z, qs, prec);

        z++;
    } while (acb_dirichlet_conrey_next(cn, G) >= 0);

    acb_dirichlet_dft_conrey(res, zeta, G, prec);

    acb_dirichlet_conrey_clear(cn);
    _acb_vec_clear(zeta, G->phi_q);
    acb_clear(qs);
    acb_clear(a);
}
