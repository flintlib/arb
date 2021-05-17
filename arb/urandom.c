/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_urandom(arb_t x, flint_rand_t state, slong prec, arf_rnd_t rnd)
{
    arf_urandom(arb_midref(x), state, prec, rnd);
    mag_zero(arb_radref(x));
}

