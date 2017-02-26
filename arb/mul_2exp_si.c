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
arb_mul_2exp_si(arb_t y, const arb_t x, slong e)
{
    arf_mul_2exp_si(arb_midref(y), arb_midref(x), e);
    mag_mul_2exp_si(arb_radref(y), arb_radref(x), e);
}

