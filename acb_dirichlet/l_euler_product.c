/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

#define ONE_OVER_LOG2 1.4426950408889634

void
acb_dirichlet_l_euler_product(acb_t res, const acb_t s,
    const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    arf_t left;
    slong wp, powprec, left_s;
    ulong p, p_limit;
    double p_needed_approx, powmag, logp, errmag;
    int is_real;
    acb_t t, u, v, c;
    mag_t err;

    if (!acb_is_finite(s))
    {
        acb_indeterminate(res);
        return;
    }

    arf_init(left);
    arf_set_mag(left, arb_radref(acb_realref(s)));
    arf_sub(left, arb_midref(acb_realref(s)), left, 2 * MAG_BITS, ARF_RND_FLOOR);

    /* Require re(s) >= 2. */
    if (arf_cmp_si(left, 2) < 0)
    {
        acb_indeterminate(res);
        arf_clear(left);
        return;
    }

    is_real = acb_is_real(s) && dirichlet_char_is_real(G, chi);

    /* L(s) ~= 1. */
    if (arf_cmp_si(left, prec) > 0)
    {
        acb_one(res);
        mag_hurwitz_zeta_uiui(arb_radref(acb_realref(res)), prec, 2);
        if (!is_real)
            mag_set(arb_radref(acb_imagref(res)), arb_radref(acb_realref(res)));
        acb_inv(res, res, prec);
        arf_clear(left);
        return;
    }

    left_s = arf_get_si(left, ARF_RND_FLOOR);

    /* Heuristic. */
    wp = prec + FLINT_BIT_COUNT(prec) + (prec / left_s) + 4;

    /* Don't work too hard if a small s was passed as input. */
    p_limit = 100 + prec * sqrt(prec);

    /* Truncating at p ~= 2^(prec/s) gives an error of 2^-prec */
    if (((double) prec) / left_s > 50.0)
        p_needed_approx = pow(2.0, 50.0);
    else
        p_needed_approx = pow(2.0, ((double) prec) / left_s);

    p_needed_approx = FLINT_MIN(p_limit, p_needed_approx);
    /* todo: use this to determine whether to precompute chi */

    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(c);

    acb_one(v);

    for (p = 2; p < p_limit; p = n_nextprime(p, 1))
    {
        /* p^s */
        logp = log(p);
        powmag = left_s * logp * ONE_OVER_LOG2;

        /* zeta(s,p) ~= 1/p^s + 1/((s-1) p^(s-1)) */
        errmag = (log(left_s - 1.0) + (left_s - 1.0) * logp) * ONE_OVER_LOG2;
        errmag = FLINT_MIN(powmag, errmag);

        if (errmag > prec + 2)
            break;

        powprec = FLINT_MAX(wp - powmag, 8);

        acb_dirichlet_chi(c, G, chi, p, powprec);

        if (!acb_is_zero(c))
        {
            acb_set_ui(t, p);
            acb_pow(t, t, s, powprec);
            acb_set_round(u, v, powprec);
            acb_div(t, u, t, powprec);
            acb_mul(t, t, c, powprec);
            acb_sub(v, v, t, wp);
        }
    }

    mag_init(err);
    mag_hurwitz_zeta_uiui(err, left_s, p);
    if (is_real)
        arb_add_error_mag(acb_realref(v), err);
    else
        acb_add_error_mag(v, err);
    mag_clear(err);

    acb_inv(res, v, prec);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(c);
    arf_clear(left);
}

