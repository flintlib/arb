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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_get_fmpz_mid_rad_10exp(fmpz_t mid, fmpz_t rad, fmpz_t exp, const arb_t x, long n)
{
    fmpz_t e, m;
    arb_t t, u;
    arf_t r;
    long prec;
    int roundmid, roundrad;

    if (!arb_is_finite(x) || arb_is_zero(x))
    {
        fmpz_zero(mid);
        fmpz_zero(rad);
        fmpz_zero(exp);
        return;
    }

    /*
      We compute m such that x * 10^m ~= 10^(n+5).
      If x = 2^e then m = (n+5) - e*log(2)/log(10).
    */
    fmpz_init(e);
    fmpz_init(m);
    arb_init(t);
    arb_init(u);
    arf_init(r);

    if (arf_cmpabs_mag(arb_midref(x), arb_radref(x)) > 0)
        fmpz_set(e, ARF_EXPREF(arb_midref(x)));
    else
        fmpz_set(e, ARF_EXPREF(arb_radref(x)));

    prec = fmpz_bits(e) + 15;

    arb_const_log2(t, prec);
    arb_const_log10(u, prec);
    arb_div(t, t, u, prec);
    arb_mul_fmpz(t, t, e, prec);
    arb_neg(t, t);
    arb_add_ui(t, t, n + 5, prec);

    arf_get_fmpz(m, arb_midref(t), ARF_RND_FLOOR);

    arb_set_ui(t, 10);

    fmpz_neg(exp, m);

    prec = n * 3.32192809488736 + 30;

    if (fmpz_sgn(m) >= 0)
    {
        arb_pow_fmpz_binexp(t, t, m, prec + 2 * fmpz_bits(m));
        arb_mul(t, x, t, prec);
    }
    else
    {
        fmpz_neg(m, m);
        arb_pow_fmpz_binexp(t, t, m, prec + 2 * fmpz_bits(m));
        arb_div(t, x, t, prec);
    }

    roundmid = arf_get_fmpz_fixed_si(mid, arb_midref(t), 0);

    arf_set_mag(r, arb_radref(t));
    roundrad = arf_get_fmpz_fixed_si(rad, r, 0);

    fmpz_add_ui(rad, rad, roundmid + roundrad);

    fmpz_clear(e);
    fmpz_clear(m);
    arb_clear(t);
    arb_clear(u);
    arf_clear(r);
}

