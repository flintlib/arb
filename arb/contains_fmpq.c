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
arb_contains_fmpq(const arb_t x, const fmpq_t q)
{
    fmpz_t a, b;
    long exp = *arb_expref(x);
    int result;

    /* fixme */
    if (COEFF_IS_MPZ(exp))
    {
        printf("error: arb_contains_fmpq: too large exponent\n");
        abort();
    }

    fmpz_init(a);
    fmpz_init(b);

    /* compare with left endpoint */
    fmpz_sub(a, arb_midref(x), arb_radref(x));
    fmpz_mul(a, a, fmpq_denref(q));
    if (exp >= 0)
    {
        fmpz_mul_2exp(a, a, exp);
        fmpz_set(b, fmpq_numref(q));
    }
    else
        fmpz_mul_2exp(b, fmpq_numref(q), -exp);

    result = (fmpz_cmp(a, b) <= 0);

    /* compare with right endpoint */
    if (result != 0)
    {
        fmpz_add(a, arb_midref(x), arb_radref(x));
        fmpz_mul(a, a, fmpq_denref(q));
        if (exp >= 0)
            fmpz_mul_2exp(a, a, exp);
        result = (fmpz_cmp(a, b) >= 0);
    }

    fmpz_clear(a);
    fmpz_clear(b);

    return result;
}
