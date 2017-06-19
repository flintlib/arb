/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_hurwitz(acb_t res, const acb_t s, const acb_t a, slong prec)
{
    if (acb_is_one(a))
    {
        acb_dirichlet_zeta(res, s, prec);
        return;
    }

    if (acb_is_zero(s))
    {
        acb_mul_2exp_si(res, a, 1);
        acb_sub_ui(res, res, 1, prec);
        acb_neg(res, res);
        acb_mul_2exp_si(res, res, -1);
        return;
    }

    _acb_poly_zeta_cpx_series(res, s, a, 0, 1, prec);
}

