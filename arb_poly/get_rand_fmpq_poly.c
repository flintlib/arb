/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

void
arb_poly_get_rand_fmpq_poly(fmpq_poly_t q, flint_rand_t state, const arb_poly_t x)
{
    long i, exp = *arb_expref(x);
    fmpz * den;

    if (COEFF_IS_MPZ(exp))
    {
        printf("error: arb_get_rand_fmpq: too large exponent\n");
        abort();
    }

    fmpq_poly_fit_length(q, x->length);
    den = fmpq_poly_denref(q);

    /* pick a denominator */
    fmpz_randbits(den, state, 1 + n_randint(state, arb_poly_prec(x)));
    fmpz_abs(den, den);
    if (fmpz_is_zero(den))
        fmpz_set_ui(den, 1UL);

    for (i = 0; i < x->length; i++)
        _arb_get_rand_fmpq(fmpq_poly_numref(q) + i, state, den,
            arb_poly_coeffs(x) + i, arb_poly_radref(x), exp);

    if (exp < 0)
    {
        fmpz_mul_2exp(den, den, -exp);
    }

    _fmpq_poly_set_length(q, x->length);
    _fmpq_poly_normalise(q);
    fmpq_poly_canonicalise(q);
}
