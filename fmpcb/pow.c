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
fmpcb_pow_fmpz_binexp(fmpcb_t y, const fmpcb_t b, const fmpz_t e, long prec)
{
    long i, wp, bits;

    if (-2L <= *e && *e <= 2L)
    {
        if (*e == 0L)
            fmpcb_one(y);
        else if (*e == 1L)
            fmpcb_set_round(y, b, prec);
        else if (*e == -1L)
            fmpcb_inv(y, b, prec);
        else if (*e == 2L)
            fmpcb_mul(y, b, b, prec);
        else
        {
            fmpcb_inv(y, b, prec);
            fmpcb_mul(y, y, y, prec);
        }
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        fmpcb_pow_fmpz_binexp(y, b, f, prec + 2);
        fmpcb_inv(y, y, prec);
        fmpz_clear(f);
        return;
    }

    if (y == b)
    {
        fmpcb_t t;
        fmpcb_init(t);
        fmpcb_set(t, b);
        fmpcb_pow_fmpz_binexp(y, t, e, prec);
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
fmpcb_pow_fmpz(fmpcb_t y, const fmpcb_t b, const fmpz_t e, long prec)
{
    fmpcb_pow_fmpz_binexp(y, b, e, prec);
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
_fmpcb_pow_exp(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)
{
    fmpcb_t t;
    fmpcb_init(t);
    fmpcb_log(t, x, prec);
    fmpcb_mul(t, t, y, prec);
    fmpcb_exp(z, t, prec);
    fmpcb_clear(t);
}

void
_fmpcb_pow_fmprb_exp(fmpcb_t z, const fmpcb_t x, const fmprb_t y, long prec)
{
    fmpcb_t t;
    fmpcb_init(t);
    fmpcb_log(t, x, prec);
    fmpcb_mul_fmprb(t, t, y, prec);
    fmpcb_exp(z, t, prec);
    fmpcb_clear(t);
}

void
fmpcb_pow_fmprb(fmpcb_t z, const fmpcb_t x, const fmprb_t y, long prec)
{
    if (fmprb_is_exact(y))
    {
        if (fmpr_is_zero(fmprb_midref(y)))
        {
            fmpcb_one(z);
            return;
        }

        if (!fmpr_is_special(fmprb_midref(y)))
        {
            const fmpz * exp_exp = fmpr_expref(fmprb_midref(y));

            /* smallish integer powers and square roots */
            if (!COEFF_IS_MPZ(*exp_exp) && (*exp_exp >= -1L))
            {
                long exp_bits;

                exp_bits = *exp_exp + fmpz_bits(fmpr_manref(fmprb_midref(y)));

                if (exp_bits < 64)
                {
                    fmpz_t e;
                    fmpz_init(e);

                    if (*exp_exp == -1L)
                    {
                        fmpz_set(e, fmpr_manref(fmprb_midref(y)));
                        fmpcb_sqrt(z, x, prec + exp_bits);
                        fmpcb_pow_fmpz(z, z, e, prec);
                    }
                    else
                    {
                        fmpz_mul_2exp(e, fmpr_manref(fmprb_midref(y)), *exp_exp);
                        fmpcb_pow_fmpz(z, x, e, prec);
                    }

                    fmpz_clear(e);
                    return;
                }
            }
        }
    }

    _fmpcb_pow_fmprb_exp(z, x, y, prec);
}

void
fmpcb_pow(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)
{
    if (fmprb_is_zero(fmpcb_imagref(y)))
    {
        fmpcb_pow_fmprb(z, x, fmpcb_realref(y), prec);
    }
    else
    {
        _fmpcb_pow_exp(z, x, y, prec);
    }
}

