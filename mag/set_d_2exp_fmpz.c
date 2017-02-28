/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_set_d_2exp_fmpz(mag_t z, double c, const fmpz_t exp)
{
    if (c == 0.0)
    {
        mag_zero(z);
    }
    else if (c > 1e300 || c < 0.0) /* not implemented */
    {
        flint_printf("mag_set_d_2exp_fmpz\n");
        flint_abort();
    }
    else
    {
        int cexp, fix;
        mp_limb_t man;

        c = frexp(c, &cexp);

        man = (mp_limb_t)(c * (double)(LIMB_ONE << MAG_BITS)) + 1;

        fix = man >> (MAG_BITS);
        man = (man >> fix) + fix;  /* XXX: need +fix? */
        MAG_MAN(z) = man;
        _fmpz_add_fast(MAG_EXPREF(z), exp, cexp + fix);
    }
}
