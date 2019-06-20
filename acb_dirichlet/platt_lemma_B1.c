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
_arb_ui_pow_arb(arb_t res, ulong n, const arb_t x, slong prec)
{
    arb_t a;
    arb_init(a);
    arb_set_ui(a, n);
    arb_pow(res, a, x, prec);
    arb_clear(a);
}

void
acb_dirichlet_platt_lemma_B1(arb_t out,
        slong sigma, const arb_t t0, const arb_t h, slong J, slong prec)
{
    arb_t pi, C;
    arb_t x1, x2, x3;

    if (sigma % 2 == 0 || sigma < 3)
    {
        arb_zero_pm_inf(out);
        return;
    }

    arb_init(pi);
    arb_init(C);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);

    arb_const_pi(pi, prec);
    acb_dirichlet_platt_c_bound(C, sigma, t0, h, 0, prec);

    arb_set_si(x1, 2*sigma - 1);
    arb_div(x1, x1, h, prec);
    arb_sqr(x1, x1, prec);
    arb_mul_2exp_si(x1, x1, -3);
    arb_exp(x1, x1, prec);

    arb_set_si(x2, 1 - 2*sigma);
    arb_mul_2exp_si(x2, x2, -2);
    arb_pow(x2, pi, x2, prec);

    arb_set_si(x3, 1 - sigma);
    _arb_ui_pow_arb(x3, (ulong) J, x3, prec);
    arb_div_si(x3, x3, sigma - 1, prec);

    arb_mul(out, x1, x2, prec);
    arb_mul(out, out, x3, prec);
    arb_mul(out, out, C, prec);

    arb_clear(pi);
    arb_clear(C);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
}
