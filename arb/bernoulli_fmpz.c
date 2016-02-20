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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_bernoulli_fmpz(arb_t res, const fmpz_t n, slong prec)
{
    if (fmpz_cmp_ui(n, UWORD_MAX) <= 0)
    {
        if (fmpz_sgn(n) >= 0)
            arb_bernoulli_ui(res, fmpz_get_ui(n), prec);
        else
            arb_zero(res);
    }
    else if (fmpz_is_odd(n))
    {
        arb_zero(res);
    }
    else
    {
        arb_t t;
        slong wp;

        arb_init(t);
        wp = prec + 2 * fmpz_bits(n);

        /* zeta(n) ~= 1 */
        arf_one(arb_midref(res));
        mag_one(arb_radref(res));
        mag_mul_2exp_si(arb_radref(res), arb_radref(res), WORD_MIN);

        /* |B_n| = 2 * Gamma(n)! / (2*pi)^n * zeta(n) */
        arb_gamma_fmpz(t, n, wp);
        arb_mul_fmpz(t, t, n, wp);
        arb_mul(res, res, t, wp);

        arb_const_pi(t, wp);
        arb_mul_2exp_si(t, t, 1);
        arb_pow_fmpz(t, t, n, wp);

        arb_div(res, res, t, prec);
        arb_mul_2exp_si(res, res, 1);

        if (fmpz_fdiv_ui(n, 4) == 0)
            arb_neg(res, res);

        arb_clear(t);
    }
}

