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
acb_dirichlet_l_hurwitz(acb_t res, const acb_t s,
    const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{
    ulong chin;
    acb_t z, t, u, a;
    acb_ptr xz;
    acb_dirichlet_conrey_t cn;

    acb_dirichlet_conrey_init(cn, G);
    acb_init(z);
    acb_init(t);
    acb_init(u);
    acb_init(a);

    acb_dirichlet_conrey_one(cn, G);
    acb_zero(t);

    acb_dirichlet_nth_root(z, chi->order.n, prec);
    xz = _acb_vec_init(chi->order.n);
    _acb_vec_set_powers(xz, z, chi->order.n, prec);

    do {
        chin = acb_dirichlet_ui_chi_conrey(G, chi, cn);

        acb_set_ui(a, cn->n);
        acb_div_ui(a, a, G->q, prec);
        acb_hurwitz_zeta(u, s, a, prec);

        acb_addmul(t, xz + chin, u, prec);

    } while (acb_dirichlet_conrey_next(cn, G) >= 0);

    acb_set_ui(u, G->q);
    acb_neg(a, s);
    acb_pow(u, u, a, prec);
    acb_mul(res, t, u, prec);

    acb_dirichlet_conrey_clear(cn);

    _acb_vec_clear(xz, chi->order.n);
    acb_clear(z);
    acb_clear(t);
    acb_clear(u);
    acb_clear(a);
}
