/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_eta(acb_t res, const acb_t s, slong prec)
{
    if (acb_is_one(s))
    {
        arb_const_log2(acb_realref(res), prec);
        arb_zero(acb_imagref(res));
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_one(t);
        acb_mul_2exp_si(t, t, -1);
        acb_pow(t, t, s, prec);
        acb_mul_2exp_si(t, t, 1);
        acb_sub_ui(t, t, 1, prec);
        acb_neg(t, t);
        acb_zeta(res, s, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
    }
}
