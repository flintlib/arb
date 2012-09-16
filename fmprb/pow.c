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

void
fmprb_pow_fmpz(fmprb_t y, const fmprb_t b, const fmpz_t e, long prec)
{
    long i, wp, bits;

    if (fmpz_is_zero(e))
    {
        fmprb_set_ui(y, 1UL);
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        fmprb_pow_fmpz(y, b, f, prec + 2);
        fmprb_ui_div(y, 1UL, y, prec);
        fmpz_clear(f);
    }

    if (y == b)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_set(t, b);
        fmprb_pow_fmpz(y, t, e, prec);
        fmprb_clear(t);
        return;
    }

    fmprb_set(y, b);

    bits = fmpz_bits(e);
    wp = FMPR_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        fmprb_mul(y, y, y, wp);
        if (fmpz_tstbit(e, i))
            fmprb_mul(y, y, b, wp);
    }
}

void
fmprb_pow_ui(fmprb_t y, const fmprb_t b, ulong e, long prec)
{
    fmpz_t f;
    fmpz_init_set_ui(f, e);
    fmprb_pow_fmpz(y, b, f, prec);
    fmpz_clear(f);
}

void
fmprb_ui_pow_ui(fmprb_t y, ulong b, ulong e, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, b);
    fmprb_pow_ui(y, t, e, prec);
    fmprb_clear(t);
}

void
fmprb_si_pow_ui(fmprb_t y, long b, ulong e, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_si(t, b);
    fmprb_pow_ui(y, t, e, prec);
    fmprb_clear(t);
}

void
fmprb_pow_fmpq(fmprb_t y, const fmprb_t x, const fmpq_t a, long prec)
{
    if (fmpz_is_one(fmpq_denref(a)))
    {
        fmprb_pow_fmpz(y, x, fmpq_numref(a), prec);
    }
    /* TODO: generalize this to a = p/q for any small p, q */
    else if (fmpz_is_one(fmpq_numref(a)) && fmpz_cmp_ui(fmpq_denref(a), 2) == 0)
    {
        fmprb_sqrt(y, x, prec);
    }
    else
    {
        long wp;

        wp = prec + 10;

        fmprb_log(y, x, wp);
        fmprb_mul_fmpz(y, y, fmpq_numref(a), wp);
        fmprb_div_fmpz(y, y, fmpq_denref(a), wp);
        fmprb_exp(y, y, prec);
    }
}
