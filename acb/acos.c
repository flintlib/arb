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
acb_acos(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_one(z))
    {
        acb_zero(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        if (arb_is_zero(acb_imagref(z)))
        {
            arb_t one;
            arb_init(one);
            arb_one(one);
            if (arb_gt(acb_realref(z), one))
            {
                acb_asin(res, z, prec);
                acb_neg(res, res);
                arb_zero(acb_realref(res));
                arb_clear(one);
                arb_clear(t);
                return;
            }
            arb_clear(one);
        }
        acb_asin(res, z, prec);
        acb_const_pi(t, prec);
        acb_mul_2exp_si(t, t, -1);
        acb_sub(res, t, res, prec);
        acb_clear(t);
    }
}

