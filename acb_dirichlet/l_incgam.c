/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "acb_hypgeom.h"

ulong
_acb_dirichlet_l_incgam_length(const acb_t s, ulong q, slong prec)
{
#if 0
    arf_t a;
    arf_init(a);
    arb_get_lbound_arf(a, s, 53);
    a = arf_get_d(at, ARF_RND_DOWN);
#else
    /* as a first approximation for Re(s) <= 1 */
    return acb_dirichlet_theta_length_d(q, 1., prec); 
#endif
}

void
acb_dirichlet_l_incgam(acb_t res, const acb_t s,
    const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{
    slong k, len;
    acb_ptr a;
    acb_t s1, s2, args, arg1s, eps;
    acb_t g1, g2, pq, pqk2;

    len = _acb_dirichlet_l_incgam_length(s, chi->q, prec) + 2;
    flint_printf("incgam length = %ld\n",len);

    prec += n_clog(len, 2);

    a = _acb_vec_init(len);
    acb_dirichlet_chi_vec(a, G, chi, len, prec);

    acb_init(args);
    acb_init(arg1s);
    acb_init(s1);
    acb_init(s2);
    acb_init(eps);
    acb_init(g1);
    acb_init(g2);
    acb_init(pq);
    acb_init(pqk2);

    acb_set(args, s);
    acb_one(arg1s);
    acb_sub(arg1s, arg1s, s, prec);

    if (chi->parity)
    {
        acb_add_ui(args, args, 1, prec);
        acb_add_ui(arg1s, arg1s, 1, prec);

        for (k = 2; k < len; k++)
            acb_mul_si(a + k, a + k, k, prec);

    }

    acb_mul_2exp_si(args, args, -1);
    acb_mul_2exp_si(arg1s, arg1s, -1);
    acb_conj(arg1s, arg1s);

    arb_const_pi(acb_realref(pq), prec);
    acb_div_ui(pq, pq, G->q, prec);

    acb_dirichlet_root_number(eps, G, chi, prec);

    for (k = 1; k < len; k++)
    {
        if (acb_is_zero(a + k))
            continue;
        acb_mul_ui(pqk2, pq, k*k, prec);
        /* FIXME: accuracy can be very bad */
        acb_hypgeom_gamma_upper(g1, args, pqk2, 2, prec);
        acb_hypgeom_gamma_upper(g2, arg1s, pqk2, 2, prec);
        acb_addmul(s1, a + k, g1, prec);
        acb_addmul(s2, a + k, g2, prec);
    }

    acb_conj(s2, s2);
    acb_mul(s2, eps, s2, prec);
    acb_add(res, s1, s2, prec);

    acb_pow(g1, pq, args, prec);
    acb_mul(res, res, g1, prec);

    acb_gamma(g1, args, prec);
    acb_div(res, res, g1, prec);

    acb_clear(args);
    acb_clear(arg1s);
    acb_clear(s1);
    acb_clear(s2);
    acb_clear(eps);
    acb_clear(g1);
    acb_clear(g2);
    acb_clear(pq);
    acb_clear(pqk2);
    _acb_vec_clear(a, len);
}
