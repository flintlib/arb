/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_clear(arf_t x)
{
    ARF_DEMOTE(x);

    /* fmpz_clear(), but avoids a redundant store */
    if (COEFF_IS_MPZ(ARF_EXP(x)))
        _fmpz_clear_mpz(ARF_EXP(x));
}

