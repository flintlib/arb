/*
    Copyright (C) 2018 D.H.J Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void acb_riemann_xi(acb_t res, const acb_t s, slong prec)
{
    acb_t pi, z1, z2, z3, z4;

    acb_init(pi);
    acb_const_pi(pi, prec);

    /* s*(s-1)/2 */
    acb_init(z1);
    acb_sub_ui(z1, s, 1, prec);
    acb_mul(z1, z1, s, prec);
    acb_mul_2exp_si(z1, z1, -1);

    /* pi^(-s/2) */
    acb_init(z2);
    acb_mul_2exp_si(z2, s, -1);
    acb_neg(z2, z2);
    acb_pow(z2, pi, z2, prec);

    /* gamma(s/2) */
    acb_init(z3);
    acb_mul_2exp_si(z3, s, -1);
    acb_gamma(z3, z3, prec);

    /* zeta(s) */
    acb_init(z4);
    acb_zeta(z4, s, prec);

    acb_mul(res, z1, z2, prec);
    acb_mul(res, res, z3, prec);
    acb_mul(res, res, z4, prec);

    acb_clear(pi);
    acb_clear(z1);
    acb_clear(z2);
    acb_clear(z3);
    acb_clear(z4);
}
