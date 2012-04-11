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
#include "ufloat.h"

/* TODO: check for overflow (both mid and rad) */


/*
exp((a+b)*2^r) - exp(a*2^r) = exp(a*2^r) * (exp(b*2^r)-1)
*/
void
arb_exp_error(ufloat_t err, const fmpz_t mid, const fmpz_t rad, const fmpz_t exp)
{
    mpfr_t a, b;

    mpfr_init2(a, 32);
    mpfr_init2(b, 32);

    _arb_get_mpfr(a, mid, exp, MPFR_RNDU);
    _arb_get_mpfr(b, rad, exp, MPFR_RNDU);

    mpfr_exp(a, a, MPFR_RNDU);
    mpfr_expm1(b, b, MPFR_RNDU);
    mpfr_mul(a, a, b, MPFR_RNDU);

    ufloat_set_mpfr(err, a);

    mpfr_clear(a);
    mpfr_clear(b);
}

void
arb_exp(arb_t y, const arb_t x)
{
    long prec;
    mpfr_t t, u;
    int input_approx;
    ufloat_t err;

    prec = FLINT_MAX(2, arb_prec(y));

    mpfr_init2(t, 2 + fmpz_bits(arb_midref(x)));
    mpfr_init2(u, prec);
    arb_get_mpfr(t, x, MPFR_RNDN);  /* exact */

    input_approx = !fmpz_is_zero(arb_radref(x));
    if (input_approx)
        arb_exp_error(err, arb_midref(x), arb_radref(x), arb_expref(x));

    mpfr_exp(u, t, MPFR_RNDN);
    arb_set_mpfr(y, u, 1);

    if (input_approx)
        _arb_rad_add_ufloat(y, err);

    mpfr_clear(t);
    mpfr_clear(u);
}
