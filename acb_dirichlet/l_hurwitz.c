/*
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "acb_poly.h"

/* todo: should document or fix that it doesn't allow aliasing */
void
acb_dirichlet_l_hurwitz(acb_t res, const acb_t s,
    const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{
    ulong chin;
    acb_t t, u, a;
    acb_ptr z;
    acb_dirichlet_conrey_t cn;
    int deflate;

    /* remove pole in Hurwitz zeta at s = 1 */
    deflate = 0;
    if (acb_is_one(s))
    {
        /* character is principal */
        if (chi->x->n == 1)
        {
            acb_indeterminate(res);
            return;
        }
        deflate = 1;
    }

    acb_dirichlet_conrey_init(cn, G);
    acb_init(t);
    acb_init(u);
    acb_init(a);

    acb_dirichlet_conrey_one(cn, G);
    acb_zero(t);

    prec += n_clog(G->phi_q, 2);

    z = _acb_vec_init(chi->order.n);
    acb_dirichlet_vec_nth_roots(z, chi->order.n, prec);

    do {
        chin = acb_dirichlet_ui_chi_conrey(G, chi, cn);

        acb_set_ui(a, cn->n);
        acb_div_ui(a, a, G->q, prec);

        if (deflate == 0)
            acb_hurwitz_zeta(u, s, a, prec);
        else
            _acb_poly_zeta_cpx_series(u, s, a, 1, 1, prec);

        acb_addmul(t, z + chin, u, prec);

    } while (acb_dirichlet_conrey_next(cn, G) >= 0);

    acb_set_ui(u, G->q);
    acb_neg(a, s);
    acb_pow(u, u, a, prec);
    acb_mul(res, t, u, prec);

    acb_dirichlet_conrey_clear(cn);

    _acb_vec_clear(z, chi->order.n);
    acb_clear(t);
    acb_clear(u);
    acb_clear(a);
}

