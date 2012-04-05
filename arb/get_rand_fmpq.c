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
arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x)
{
    fmpz_t a, b;
    long exp = *arb_expref(x);

    if (COEFF_IS_MPZ(exp))
    {
        printf("error: arb_get_rand_fmpq: too large exponent\n");
        abort();
    }

    /* there is only one rational */
    if (fmpz_is_zero(arb_radref(x)))
    {
        fmpz_set(fmpq_numref(q), arb_midref(x));
        fmpz_set_ui(fmpq_denref(q), 1UL);
        if (exp >= 0)
            fmpz_mul_2exp(fmpq_numref(q), fmpq_numref(q), exp);
        else
            fmpz_mul_2exp(fmpq_denref(q), fmpq_denref(q), -exp);
        fmpq_canonicalise(q);
        return;
    }

    fmpz_init(a);
    fmpz_init(b);

    fmpz_sub(a, arb_midref(x), arb_radref(x));
    fmpz_add(b, arb_midref(x), arb_radref(x));

    /* pick a denominator */
    fmpz_randbits(fmpq_denref(q), state, 1 + n_randint(state, arb_prec(x)));
    fmpz_abs(fmpq_denref(q), fmpq_denref(q));
    if (fmpz_is_zero(fmpq_denref(q)))
        fmpz_set_ui(fmpq_denref(q), 1UL);

    if (exp >= 0)
    {
        fmpz_mul_2exp(a, a, exp);
        fmpz_mul_2exp(b, b, exp);
    }

    /* generate random integer in [a*den, b*den] */
    fmpz_mul(a, a, fmpq_denref(q));
    fmpz_mul(b, b, fmpq_denref(q));
    fmpz_add_ui(b, b, 1UL);
    fmpz_sub(b, b, a);
    fmpz_randtest_mod(fmpq_numref(q), state, b);
    fmpz_add(fmpq_numref(q), fmpq_numref(q), a);

    if (exp < 0)
        fmpz_mul_2exp(fmpq_denref(q), fmpq_denref(q), -exp);

    fmpq_canonicalise(q);

    fmpz_clear(a);
    fmpz_clear(b);
}
