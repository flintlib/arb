/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
_acb_vec_nth_roots(acb_ptr res, slong n, slong prec)
{
    slong k, len, wp;
    acb_t t;

    if (n <= 0)
        return;

    if (n % 4 == 0)
        len = n / 8 + 1;
    else if (n % 2 == 0)
        len = n / 4 + 1;
    else
        len = n / 2 + 1;

    wp = prec + 6 + 2 * FLINT_BIT_COUNT(len);

    acb_init(t);
    acb_nth_root(t, n, wp);
    _acb_vec_set_powers(res, t, len, wp);
    acb_clear(t);
    _acb_vec_set_round(res, res, len, prec);

    if (n % 4 == 0)
    {
        for (k = n / 8 + 1; k <= n / 4; k++)
        {
            arb_set(acb_realref(res + k), acb_imagref(res + n / 4 - k));
            arb_set(acb_imagref(res + k), acb_realref(res + n / 4 - k));
        }

        for (k = n / 4 + 1; k <= n / 2; k++)
            acb_mul_onei(res + k, res + k - n / 4);
    }
    else if (n % 2 == 0)
    {
        for (k = n / 4 + 1; k <= n / 2; k++)
        {
            acb_set(res + k, res + n / 2 - k);
            arb_neg(acb_realref(res + k), acb_realref(res + k));
        }
    }

    for (k = n / 2 + 1; k < n; k++)
        acb_conj(res + k, res + n - k);
}

