/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static void
_arb_div_ui_ui(arb_t res, ulong a, ulong b, slong prec)
{
    arb_set_ui(res, a);
    arb_div_ui(res, res, b, prec);
}

/*
 * This function implements Theorem 1 and both parts of (1.2) from
 * "An improved upper bound for the argument of the Riemann
 * zeta-function on the critical line II"
 * by Timothy Trudgian.
 */
static void
_backlund_s_trudgian_bound(mag_t res, const arb_t t, slong prec)
{
    if (!arb_is_nonnegative(t))
    {
        mag_inf(res);
    }
    else
    {
        arf_t u;
        arf_init(u);
        arb_get_ubound_arf(u, t, prec);
        if (arf_cmp_si(u, 280) < 0)
        {
            mag_one(res);
        }
        else if (arf_cmp_si(u, 6800000) < 0)
        {
            mag_set_ui(res, 2);
        }
        else
        {
            /* |S(u)| <= 0.111*log(u) + 0.275*log(log(u)) + 2.45 */
            arb_t c, r, logu;
            arb_init(c);
            arb_init(r);
            arb_init(logu);
            arb_log_arf(logu, u, prec);
            _arb_div_ui_ui(c, 275, 1000, prec);
            arb_log(r, logu, prec);
            arb_mul(r, r, c, prec);
            _arb_div_ui_ui(c, 111, 1000, prec);
            arb_addmul(r, c, logu, prec);
            _arb_div_ui_ui(c, 245, 100, prec);
            arb_add(r, r, c, prec);
            arb_get_mag(res, r);
            arb_clear(c);
            arb_clear(r);
            arb_clear(logu);
        }
        arf_clear(u);
    }
}

void
acb_dirichlet_backlund_s(arb_t res, const arb_t t, slong prec)
{
    /* TODO: Use Turing's method for greater accuracy. */
    arb_zero(res);
    _backlund_s_trudgian_bound(arb_radref(res), t, prec);
}
