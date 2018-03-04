/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_acosh(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_one(z))
    {
        acb_zero(res);
    }
    else
    {
        acb_t t, u;
        acb_init(t);
        acb_init(u);

        acb_add_ui(t, z, 1, prec);
        acb_sub_ui(u, z, 1, prec);
        acb_sqrt(t, t, prec);
        acb_sqrt(u, u, prec);
        acb_mul(t, t, u, prec);
        acb_add(t, t, z, prec);
        acb_log(res, t, prec);

        if (arb_is_zero(acb_imagref(z)))
        {
            arb_t one, abs;
            arb_init(one);
            arb_init(abs);
            arb_one(one);
            arb_abs(abs, acb_realref(z));
            if (arb_lt(abs, one))
            {
                arb_zero(acb_realref(res));
            }
            arb_clear(one);
            arb_clear(abs);
        }

        acb_clear(t);
        acb_clear(u);
    }
}

