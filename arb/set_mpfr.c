/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb.h"

/* set to value of mpfr with given ulp error */
void
arb_set_mpfr(arb_t y, const mpfr_t x, ulong error)
{
    if (mpfr_zero_p(x))
    {
        arb_zero(y);
        return;
    }
    else
    {
        mpz_t t;
        long exp;

        mpz_init(t);

        exp = mpfr_get_z_2exp(t, x);

        fmpz_set_mpz(arb_midref(y), t);
        fmpz_set_ui(arb_radref(y), error);
        fmpz_set_si(arb_expref(y), exp);

        mpz_clear(t);
    }
}
