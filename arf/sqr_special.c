/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_sqr_special(arf_t res, const arf_t x)
{
    if (arf_is_zero(x))
        arf_zero(res);
    else if (arf_is_nan(x))
        arf_nan(res);
    else
        arf_pos_inf(res);
}
