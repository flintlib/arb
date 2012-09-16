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

#include "fmpr.h"

/* TODO: correct rounding with negative b */

void
fmpr_pow_sloppy_fmpz(fmpr_t y, const fmpr_t b, const fmpz_t e,
    long prec, fmpr_rnd_t rnd)
{
    long i, wp, bits;

    if (fmpz_is_zero(e))
    {
        fmpr_set_ui(y, 1UL);
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        fmpr_pow_sloppy_fmpz(y, b, f, prec + 2,
            (rnd == FMPR_RND_FLOOR || rnd == FMPR_RND_DOWN)
            ? FMPR_RND_UP : FMPR_RND_DOWN);
        fmpr_ui_div(y, 1UL, y, prec, rnd);
        fmpz_clear(f);
    }

    if (y == b)
    {
        fmpr_t t;
        fmpr_init(t);
        fmpr_set(t, b);
        fmpr_pow_sloppy_fmpz(y, t, e, prec, rnd);
        fmpr_clear(t);
        return;
    }

    fmpr_set(y, b);

    bits = fmpz_bits(e);
    wp = FMPR_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        fmpr_mul(y, y, y, wp, rnd);
        if (fmpz_tstbit(e, i))
            fmpr_mul(y, y, b, wp, rnd);
    }
}

void
fmpr_pow_sloppy_ui(fmpr_t y, const fmpr_t b, ulong e, long prec, fmpr_rnd_t rnd)
{
    fmpz_t f;
    fmpz_init_set_ui(f, e);
    fmpr_pow_sloppy_fmpz(y, b, f, prec, rnd);
    fmpz_clear(f);
}

void
fmpr_pow_sloppy_si(fmpr_t y, const fmpr_t b, long e, long prec, fmpr_rnd_t rnd)
{
    fmpz_t f;
    fmpz_init(f);
    fmpz_set_si(f, e);
    fmpr_pow_sloppy_fmpz(y, b, f, prec, rnd);
    fmpz_clear(f);
}
