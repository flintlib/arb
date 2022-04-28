/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_pow_si(arb_t y, const arb_t b, slong e, slong prec)
{
    if (e >= 0)
        arb_pow_ui_binexp(y, b, e, prec);
    else if (e == -1)
        arb_inv(y, b, prec);
    else if (e == -2)
    {
        arb_inv(y, b, prec);
        arb_sqr(y, y, prec);
    }
    else if (arb_is_exact(b))
    {
        arb_pow_ui_binexp(y, b, -e, prec + 2);
        arb_inv(y, y, prec);
    }
    else
    {
        e = -e;
        arb_inv(y, b, prec + FLINT_BIT_COUNT(e) + 2);
        arb_pow_ui_binexp(y, y, e, prec);
    }
}
