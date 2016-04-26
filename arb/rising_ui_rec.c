/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_rising_ui_rec(arb_t y, const arb_t x, ulong n, slong prec)
{
    if (prec < 512 || n < 8 || arb_bits(x) < prec / 8)
        arb_rising_ui_bs(y, x, n, prec);
    else
        arb_rising_ui_rs(y, x, n, 0, prec);
}

