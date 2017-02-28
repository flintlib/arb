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

void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s,
    const signed char * chi, int mod, int reciprocal, slong prec)
{
    slong wp, powprec;
    double logp, powmag, errmag, limit;
    arb_t t, u;
    ulong p;
    mag_t err;

    if (s <= 1)
    {
        arb_indeterminate(res);
        return;
    }

    if (prec < 2) flint_abort(); /* assert */

    /* L(s), 1/L(s) = 1 + ...  For s >= 3, zeta(s,2) < 2^(1-s). */
    if (s > (ulong) prec)
    {
        arf_one(arb_midref(res));
        mag_set_ui_2exp_si(arb_radref(res), 1, -prec);
        return;
    }

    /* L(s), 1/L(s) = 1 +/- chi(2) 2^(-s) + ...
       For s >= 2, zeta(s,3) < 2^(2-floor(3s/2)). */
    if (s > 0.7 * prec)
    {
        arb_one(res);

        if (chi[2 % mod] != 0)
        {
            arf_t t;
            arf_init(t);
            arf_set_si_2exp_si(t, chi[2 % mod], -s);
            if (reciprocal)
                arf_neg(t, t);
            arb_add_arf(res, res, t, prec);
            arf_clear(t);
        }

        arb_add_error_2exp_si(res, 2 - (3 * s) / 2);
        return;
    }

    /* Heuristic. */
    wp = prec + FLINT_BIT_COUNT(prec) + (prec / s) + 4;

    arb_init(t);
    arb_init(u);

    /* 1 - chi(2) 2^(-s) */
    arb_one(res);
    arf_set_ui_2exp_si(arb_midref(t), 1, -s);   /* -s cannot overflow */
    if (chi[2 % mod] == 1)
        arb_sub(res, res, t, wp);
    else if (chi[2 % mod] == -1)
        arb_add(res, res, t, wp);

    /* Cut at some point if this algorithm just isn't a good fit... */
    /* The C * prec / log(prec) cutoff in arb_zeta_ui implies that the limit
       should be at least prec ^ (log(2) / C) for the Riemann zeta function,
       which gives prec ^ 1.2956 here. */
    limit = 100 + prec * sqrt(prec);

    for (p = 3; p < limit; p = n_nextprime(p, 1))
    {
        if (mod == 1 || chi[p % mod] != 0)
        {
            /* p^s */
            logp = log(p);
            powmag = s * logp * ONE_OVER_LOG2;

            /* zeta(s,p) ~= 1/p^s + 1/((s-1) p^(s-1)) */
            errmag = (log(s - 1.0) + (s - 1.0) * logp) * ONE_OVER_LOG2;
            errmag = FLINT_MIN(powmag, errmag);

            if (errmag > prec + 2)
                break;

            powprec = FLINT_MAX(wp - powmag, 8);

            arb_ui_pow_ui(t, p, s, powprec);
            arb_set_round(u, res, powprec);
            arb_div(t, u, t, powprec);

            if (mod == 1 || (chi[p % mod] == 1))
                arb_sub(res, res, t, wp);
            else
                arb_add(res, res, t, wp);
        }
    }

    mag_init(err);
    mag_hurwitz_zeta_uiui(err, s, p);
    arb_add_error_mag(res, err);
    mag_clear(err);

    if (reciprocal)
        arb_set_round(res, res, prec);
    else
        arb_inv(res, res, prec);

    arb_clear(t);
    arb_clear(u);
}

