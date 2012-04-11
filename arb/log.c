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

/*
Let the input be [a-b, a+b] * 2^e. We require a > b >= 0 (otherwise the
interval contains zero or a negative number and the logarithm is not
defined). The error is largest at a-b, and we have

log(a * 2^e) - log((a-b) * 2^e) = log(1 + b/(a-b)).
*/

void
arb_log_error(ufloat_t err, const fmpz_t mid, const fmpz_t rad)
{
    ufloat_t t, u;
    fmpz_t d;

    fmpz_init(d);
    fmpz_sub(d, mid, rad);

    ufloat_set_fmpz(t, rad);
    ufloat_set_fmpz_lower(u, d);
    ufloat_div(t, t, u);
    ufloat_log1p(err, t);

    fmpz_clear(d);
}

static void
_arb_rad_add_ufloat(arb_t y, ufloat_t err)
{
    ufloat_t w;

    err->exp -= *arb_expref(y);
    ufloat_set_fmpz(w, arb_radref(y));
    ufloat_add(w, w, err);
    ufloat_get_fmpz(arb_radref(y), w);
}

void
arb_log(arb_t y, const arb_t x)
{
    long prec;
    mpfr_t t, u;
    int input_approx, value_approx;
    ufloat_t err;

    if (fmpz_sgn(arb_midref(x)) <= 0 ||
        (fmpz_cmpabs(arb_midref(x), arb_radref(x)) <= 0))
    {
        printf("arb_log: interval contains zero or negative numbers\n");
        abort();
    }

    prec = FLINT_MAX(2, arb_prec(y));

    mpfr_init2(t, 2 + fmpz_bits(arb_midref(x)));
    mpfr_init2(u, prec);
    arb_get_mpfr(t, x, MPFR_RNDN);  /* exact */

    input_approx = !fmpz_is_zero(arb_radref(x));
    if (input_approx)
        arb_log_error(err, arb_midref(x), arb_radref(x));

    value_approx = mpfr_log(u, t, MPFR_RNDN) != 0;
    arb_set_mpfr(y, u, value_approx);

    if (input_approx)
        _arb_rad_add_ufloat(y, err);

    mpfr_clear(t);
    mpfr_clear(u);
}
