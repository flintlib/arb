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
acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, slong prec)
{
    if (prec < 256 || n < 8 || acb_bits(x) < prec / 8)
        acb_rising2_ui_bs(u, v, x, n, prec);
    else
        acb_rising2_ui_rs(u, v, x, n, 0, prec);
}

