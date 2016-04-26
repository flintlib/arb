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
acb_dirichlet_l(acb_t res, const acb_t s,
    const acb_dirichlet_group_t G, ulong m, slong prec)
{
    acb_t chi, t, u, a;
    ulong k;

    acb_init(chi);
    acb_init(t);
    acb_init(u);
    acb_init(a);

    acb_zero(t);

    for (k = 1; k <= G->q; k++)
    {
        acb_dirichlet_chi(chi, G, m, k, prec);

        if (!acb_is_zero(chi))
        {
            acb_set_ui(a, k);
            acb_div_ui(a, a, G->q, prec);
            acb_hurwitz_zeta(u, s, a, prec);
            acb_addmul(t, chi, u, prec);
        }
    }

    acb_set_ui(u, G->q);
    acb_neg(a, s);
    acb_pow(u, u, a, prec);
    acb_mul(res, t, u, prec);

    acb_clear(chi);
    acb_clear(t);
    acb_clear(u);
    acb_clear(a);
}

