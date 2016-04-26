/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_inv(acb_t z, const acb_t x, slong prec)
{
#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(z)
#define d acb_imagref(z)

    if (arb_is_zero(b))
    {
        arb_inv(c, a, prec);
        arb_zero(d);
    }
    else if (arb_is_zero(a))
    {
        arb_inv(d, b, prec);
        arb_neg(d, d);
        arb_zero(c);
    }
    else
    {
        arb_t t;
        arb_init(t);

        arb_mul(t, a, a, prec);
        arb_addmul(t, b, b, prec);

        arb_div(c, a, t, prec);
        arb_div(d, b, t, prec);

        arb_neg(d, d);

        arb_clear(t);
    }

#undef a
#undef b
#undef c
#undef d
}

