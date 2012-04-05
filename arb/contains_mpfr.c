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

int
arb_contains_mpfr(const arb_t x, const mpfr_t z)
{
    fmpz_t a;
    mpfr_t t;
    long exp = *arb_expref(x);
    int result;

    /* fixme */
    if (COEFF_IS_MPZ(exp))
    {
        printf("error: arb_contains_mpfr: too large exponent\n");
        abort();
    }

    fmpz_init(a);

    /* compare with left endpoint */
    fmpz_sub(a, arb_midref(x), arb_radref(x));

    mpfr_init2(t, 2 + fmpz_bits(a));
    _arb_get_mpfr(t, a, arb_expref(x), MPFR_RNDN);  /* exact */

    result = (mpfr_cmp(t, z) <= 0);

    /* compare with right endpoint */
    if (result != 0)
    {
        fmpz_add(a, arb_midref(x), arb_radref(x));
        mpfr_set_prec(t, 2 + fmpz_bits(a));
        _arb_get_mpfr(t, a, arb_expref(x), MPFR_RNDN);  /* exact */
        result = (mpfr_cmp(t, z) >= 0);
    }

    fmpz_clear(a);
    mpfr_clear(t);
    return result;
}
