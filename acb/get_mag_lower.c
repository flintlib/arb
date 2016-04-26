/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_get_mag_lower(mag_t z, const acb_t x)
{
    arf_t t;
    arf_init(t);
    acb_get_abs_lbound_arf(t, x, MAG_BITS);
    arf_get_mag_lower(z, t);
    arf_clear(t);
}

