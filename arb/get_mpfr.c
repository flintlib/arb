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

static __inline__ void
mpfr_set_fmpz(mpfr_t c, const fmpz_t b, mpfr_rnd_t rnd)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_set_z(c, COEFF_TO_PTR(*b), rnd);
    else
        mpfr_set_si(c, *b, rnd);
}

void
_arb_get_mpfr(mpfr_t f, const fmpz_t mid, const fmpz_t exp, mpfr_rnd_t rnd)
{
    long e;

    if (COEFF_IS_MPZ(*exp))
    {
        printf("arb_get_mpfr: large exponent\n");
        abort();
    }

    e = fmpz_get_si(exp);

    mpfr_set_fmpz(f, mid, rnd);

    if (e >= 0)
        mpfr_mul_2exp(f, f, e, rnd);
    else
        mpfr_div_2exp(f, f, -e, rnd);
}

void
arb_get_mpfr(mpfr_t f, const arb_t x, mpfr_rnd_t rnd)
{
    _arb_get_mpfr(f, arb_midref(x), arb_expref(x), rnd);
}
