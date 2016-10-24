/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* todo: decompose s into midpoint+radius */
void
acb_dirichlet_zeta_rs(acb_t res, const acb_t s, slong K, slong prec)
{
    acb_t R1, R2, X, t;
    slong wp;

    if (arf_sgn(arb_midref(acb_imagref(s))) < 0)
    {
        acb_init(t);
        acb_conj(t, s);
        acb_dirichlet_zeta_rs(res, t, K, prec);
        acb_conj(res, res);
        acb_clear(t);
        return;
    }

    acb_init(R1);
    acb_init(R2);
    acb_init(X);
    acb_init(t);

    /* rs_r increases the precision internally */
    wp = prec;

    acb_dirichlet_zeta_rs_r(R1, s, K, wp);

    if (arb_is_exact(acb_realref(s)) &&
        (arf_cmp_2exp_si(arb_midref(acb_realref(s)), -1) == 0))
    {
        acb_conj(R2, R1);
    }
    else
    {
        /* conj(R(conj(1-s))) */
        arb_sub_ui(acb_realref(t), acb_realref(s), 1, 10 * wp);
        arb_neg(acb_realref(t), acb_realref(t));
        arb_set(acb_imagref(t), acb_imagref(s));
        acb_dirichlet_zeta_rs_r(R2, t, K, wp);
        acb_conj(R2, R2);
    }

    if (acb_is_finite(R1) && acb_is_finite(R2))
    {
        wp += 10 + arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(s)));
        wp = FLINT_MAX(wp, 10);

        /* X = pi^(s-1/2) gamma((1-s)/2) rgamma(s/2)
             = (2 pi)^s rgamma(s) / (2 cos(pi s / 2)) */
        acb_rgamma(X, s, wp);
        acb_const_pi(t, wp);
        acb_mul_2exp_si(t, t, 1);
        acb_pow(t, t, s, wp);
        acb_mul(X, X, t, wp);
        acb_mul_2exp_si(t, s, -1);
        acb_cos_pi(t, t, wp);
        acb_mul_2exp_si(t, t, 1);
        acb_div(X, X, t, wp);

        acb_mul(R2, R2, X, wp);
    }

    /* R1 + X * R2 */
    acb_add(res, R1, R2, prec);

    acb_clear(R1);
    acb_clear(R2);
    acb_clear(X);
    acb_clear(t);
}

