/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"

void
acb_zeta_si(acb_t z, slong s, slong prec)
{
    if (s >= 0)
    {
        arb_zeta_ui(acb_realref(z), s, prec);
    }
    else
    {
        arb_bernoulli_ui(acb_realref(z), 1-s, prec);
        arb_div_ui(acb_realref(z), acb_realref(z), 1-s, prec);
        arb_neg(acb_realref(z), acb_realref(z));
    }

    arb_zero(acb_imagref(z));
    return;
}

void
acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, slong prec)
{
    if (acb_is_one(a) && acb_is_int(s) &&
        arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), FLINT_BITS - 1) < 0)
    {
        acb_zeta_si(z, arf_get_si(arb_midref(acb_realref(s)), ARF_RND_DOWN), prec);
        return;
    }

    _acb_poly_zeta_cpx_series(z, s, a, 0, 1, prec);
}

void
acb_zeta(acb_t z, const acb_t s, slong prec)
{
    acb_dirichlet_zeta(z, s, prec);
}

