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

void
arb_log_ui(arb_t x, ulong n)
{
    mpfr_t y;
    mpz_t z;
    long exp;

    if (n == 0)
    {
        printf("arb_log_ui: log(0)\n");
        abort();
    }

    if (n == 1)
    {
        arb_zero(x);
        return;
    }

    /* power of two */
    if (n > 2 && (n & (n-1UL)) == 0UL)
    {
        arb_log_ui(x, 2);
        arb_mul_si(x, x, FLINT_BIT_COUNT(n) - 1);
        return;
    }

    mpfr_init2(y, FLINT_MAX(arb_prec(x), FLINT_BITS));

    if (n == 2)
    {
        mpfr_const_log2(y, MPFR_RNDN);
        fmpz_set_ui(arb_radref(x), 1UL);
    }
    else
    {
        mpfr_set_ui(y, n, MPFR_RNDN);  /* exact */
        mpfr_log(y, y, MPFR_RNDN);
        fmpz_set_ui(arb_radref(x), 1UL);
    }

    mpz_init(z);
    exp = mpfr_get_z_2exp(z, y);
    fmpz_set_mpz(arb_midref(x), z);
    fmpz_set_si(arb_expref(x), exp);
    mpz_clear(z);
    mpfr_clear(y);
}
