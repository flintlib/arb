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

void
acb_dirichlet_l_jet(acb_ptr res, const acb_t s,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    int deflate, slong len, slong prec)
{
    ulong order, chin, mult, phi;
    acb_t a;
    acb_ptr t, u, z;
    dirichlet_char_t cn;
    int deflate_hurwitz;

    if (len <= 0)
        return;

    if (G->q == 1)
    {
        acb_init(a);
        acb_one(a);
        _acb_poly_zeta_cpx_series(res, s, a, deflate, len, prec);
        acb_clear(a);
        return;
    }

    if (len == 1 && !(deflate && dirichlet_char_is_principal(G, chi)))
    {
        acb_dirichlet_l(res, s, G, chi, prec);
        return;
    }

    if (dirichlet_char_is_principal(G, chi))
        deflate_hurwitz = deflate;
    else
        deflate_hurwitz = acb_is_one(s);

    dirichlet_char_init(cn, G);
    t = _acb_vec_init(len);
    u = _acb_vec_init(len + 2);
    acb_init(a);

    dirichlet_char_one(cn, G);

    prec += n_clog(G->phi_q, 2);

    order = dirichlet_order_char(G, chi);
    mult = G->expo / order;
    z = _acb_vec_init(order);
    _acb_vec_nth_roots(z, order, prec);

    phi = 0;
    do
    {
        chin = dirichlet_pairing_char(G, chi, cn) / mult;
        acb_set_ui(a, cn->n);
        acb_div_ui(a, a, G->q, prec);
        _acb_poly_zeta_cpx_series(u, s, a, deflate_hurwitz, len, prec);
        _acb_vec_scalar_addmul(t, u, len, z + chin, prec);
        phi++;
    }
    while (dirichlet_char_next(cn, G) >= 0);

    if (dirichlet_char_is_principal(G, chi) && deflate)
    {
        /* res = t * q^(-(s+x)) + [phi(q) * (q^(-(s+x)) - q^-1) / ((s+x)-1)] */
        if (acb_is_one(s))
        {
            acb_set_ui(a, G->q);
            _acb_poly_acb_invpow_cpx(u, a, s, len + 1, prec);
            _acb_poly_mullow(res, t, len, u, len, len, prec);
            acb_set_ui(u, phi);
            _acb_vec_scalar_addmul(res, u + 1, len, u, prec);
        }
        else
        {
            acb_sub_ui(u, s, 1, prec);
            acb_one(u + 1);

            acb_set_ui(a, G->q);
            _acb_poly_acb_invpow_cpx(u + 2, a, s, len, prec);
            _acb_poly_mullow(res, t, len, u + 2, len, len, prec);

            acb_inv(a, a, prec);
            acb_sub(u + 2, u + 2, a, prec);
            _acb_poly_div_series(t, u + 2, len, u, 2, len, prec);

            acb_set_ui(u, phi);
            _acb_vec_scalar_addmul(res, t, len, u, prec);
        }
    }
    else
    {
        /* res = t * q^(-(s+x)) */
        acb_set_ui(a, G->q);
        _acb_poly_acb_invpow_cpx(u, a, s, len, prec);
        _acb_poly_mullow(res, t, len, u, len, len, prec);
    }

    dirichlet_char_clear(cn);

    _acb_vec_clear(z, order);
    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len + 2);
    acb_clear(a);
}

