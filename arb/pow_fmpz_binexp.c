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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_pow_fmpz_binexp(arb_t y, const arb_t b, const fmpz_t e, long prec)
{
    long i, wp, bits;

    if (-2L <= *e && *e <= 2L)
    {
        if (*e == 0L)
            arb_one(y);
        else if (*e == 1L)
            arb_set_round(y, b, prec);
        else if (*e == -1L)
            arb_inv(y, b, prec);
        else if (*e == 2L)
            arb_mul(y, b, b, prec);
        else
        {
            arb_inv(y, b, prec);
            arb_mul(y, y, y, prec);
        }
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        arb_pow_fmpz_binexp(y, b, f, prec + 2);
        arb_inv(y, y, prec);
        fmpz_clear(f);
        return;
    }

    if (y == b)
    {
        arb_t t;
        arb_init(t);
        arb_set(t, b);
        arb_pow_fmpz_binexp(y, t, e, prec);
        arb_clear(t);
        return;
    }

    arb_set(y, b);

    bits = fmpz_bits(e);
    wp = ARF_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        arb_mul(y, y, y, wp);
        if (fmpz_tstbit(e, i))
            arb_mul(y, y, b, wp);
    }
}

