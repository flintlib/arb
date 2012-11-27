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

#include "fmpcb.h"

void
fmpcb_pow_fmpz(fmpcb_t y, const fmpcb_t b, const fmpz_t e, long prec)
{
    long i, wp, bits;

    if (fmpz_is_zero(e))
    {
        fmpcb_one(y);
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        fmpcb_pow_fmpz(y, b, f, prec + 2);
        fmpcb_inv(y, y, prec);
        fmpz_clear(f);
        return;
    }

    if (y == b)
    {
        fmpcb_t t;
        fmpcb_init(t);
        fmpcb_set(t, b);
        fmpcb_pow_fmpz(y, t, e, prec);
        fmpcb_clear(t);
        return;
    }

    fmpcb_set(y, b);

    bits = fmpz_bits(e);
    wp = FMPR_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        fmpcb_mul(y, y, y, wp);
        if (fmpz_tstbit(e, i))
            fmpcb_mul(y, y, b, wp);
    }
}

void
fmpcb_pow_ui(fmpcb_t y, const fmpcb_t b, ulong e, long prec)
{
    fmpz_t f;
    fmpz_init_set_ui(f, e);
    fmpcb_pow_fmpz(y, b, f, prec);
    fmpz_clear(f);
}

void
fmpcb_pow_si(fmpcb_t y, const fmpcb_t b, long e, long prec)
{
    fmpz_t f;
    fmpz_init(f);
    fmpz_set_si(f, e);
    fmpcb_pow_fmpz(y, b, f, prec);
    fmpz_clear(f);
}

void
fmpcb_pow(fmpcb_t r, const fmpcb_t x, const fmpcb_t y, long prec)
{
    fmpcb_t t;
    fmpcb_init(t);
    fmpcb_log(t, x, prec);
    fmpcb_mul(t, t, y, prec);
    fmpcb_exp(r, t, prec);
    fmpcb_clear(t);
}

