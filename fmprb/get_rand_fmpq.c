/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

/*
Returns a random rational in x. The binary denominator of x is multiplied
by den_mult to select a base denominator to choose random integers
between.

The outcome is undefined if the midpoint or radius of x is non-finite,
or if the exponent of the midpoint or radius is so large or small
that representing the endpoints as exact rational numbers would
cause overflows.
*/

void
_fmprb_get_rand_fmpq(fmpz_t num, fmpz_t den, flint_rand_t state,
    const fmpz_t den_mult, const fmprb_t x)
{
    fmpz_t a, b, exp;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(exp);

    fmprb_get_interval_fmpz_2exp(a, b, exp, x);

    if (COEFF_IS_MPZ(*exp))
    {
        printf("exception: fmprb_get_rand_fmpq: too large exponent\n");
        abort();
    }

    if (*exp >= 0)
    {
        fmpz_mul_2exp(a, a, *exp);
        fmpz_mul_2exp(b, b, *exp);
    }

    /* generate random integer in [a*den, b*den] */
    fmpz_mul(a, a, den_mult);
    fmpz_mul(b, b, den_mult);
    fmpz_add_ui(b, b, 1UL);
    fmpz_sub(b, b, a);

    /* return one endpoint with high probability (used for stress
       testing rounding) */
    if (n_randint(state, 6) == 0)
    {
        if (n_randint(state, 2))
            fmpz_zero(num);
        else
            fmpz_sub_ui(num, b, 1UL);
    }
    else
    {
        fmpz_randtest_mod(num, state, b);
    }

    fmpz_add(num, num, a);

    fmpz_set(den, den_mult);

    if (*exp < 0)
        fmpz_mul_2exp(den, den, -(*exp));

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(exp);
}


/*
Chooses a random rational number from the interval represented by x.
A denominator is chosen by multiplying the binary denominator of x
by a random integer up to size bits.

The outcome is undefined if the midpoint or radius of x is non-finite,
or if the exponent of the midpoint or radius is so large or small
that representing the endpoints as exact rational numbers would
cause overflows.
*/

void
fmprb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const fmprb_t x, long bits)
{
    /* there is only one rational */
    if (fmprb_is_exact(x))
    {
        fmpr_get_fmpq(q, fmprb_midref(x));
        return;
    }

    /* pick a denominator */
    fmpz_randbits(fmpq_denref(q), state, n_randint(state, bits + 1));
    fmpz_abs(fmpq_denref(q), fmpq_denref(q));
    if (fmpz_is_zero(fmpq_denref(q)))
        fmpz_one(fmpq_denref(q));

    _fmprb_get_rand_fmpq(fmpq_numref(q), fmpq_denref(q), state, fmpq_denref(q), x);
    fmpq_canonicalise(q);
}
