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

#include "acb.h"

void
acb_pow_fmpz_binexp(acb_t y, const acb_t b, const fmpz_t e, long prec)
{
    long i, wp, bits;

    if (-2L <= *e && *e <= 4L)
    {
        if (*e == 0L)
        {
            acb_one(y);
        }
        else if (*e == 1L)
        {
            acb_set_round(y, b, prec);
        }
        else if (*e == -1L)
        {
            acb_inv(y, b, prec);
        }
        else if (*e == 2L)
        {
            acb_mul(y, b, b, prec);
        }
        else if (*e == 3L)
        {
            acb_cube(y, b, prec);
        }
        else if (*e == 4L)
        {
            acb_mul(y, b, b, prec);
            acb_mul(y, y, y, prec);
        }
        else
        {
            acb_inv(y, b, prec);
            acb_mul(y, y, y, prec);
        }
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        acb_pow_fmpz_binexp(y, b, f, prec + 2);
        acb_inv(y, y, prec);
        fmpz_clear(f);
        return;
    }

    if (!COEFF_IS_MPZ(*e) && ((*e) % 3 == 0))
    {
        fmpz e3 = (*e) / 3;
        acb_pow_fmpz_binexp(y, b, &e3, prec);
        acb_cube(y, y, prec);
        return;
    }

    if (y == b)
    {
        acb_t t;
        acb_init(t);
        acb_set(t, b);
        acb_pow_fmpz_binexp(y, t, e, prec);
        acb_clear(t);
        return;
    }

    acb_set(y, b);

    bits = fmpz_bits(e);
    wp = ARF_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        acb_mul(y, y, y, wp);
        if (fmpz_tstbit(e, i))
            acb_mul(y, y, b, wp);
    }
}

void
acb_pow_fmpz(acb_t y, const acb_t b, const fmpz_t e, long prec)
{
    acb_pow_fmpz_binexp(y, b, e, prec);
}

void
acb_pow_ui(acb_t y, const acb_t b, ulong e, long prec)
{
    fmpz_t f;
    fmpz_init_set_ui(f, e);
    acb_pow_fmpz(y, b, f, prec);
    fmpz_clear(f);
}

void
acb_pow_si(acb_t y, const acb_t b, long e, long prec)
{
    fmpz_t f;
    fmpz_init(f);
    fmpz_set_si(f, e);
    acb_pow_fmpz(y, b, f, prec);
    fmpz_clear(f);
}

void
_acb_pow_exp(acb_t z, const acb_t x, const acb_t y, long prec)
{
    acb_t t;
    acb_init(t);
    acb_log(t, x, prec);
    acb_mul(t, t, y, prec);
    acb_exp(z, t, prec);
    acb_clear(t);
}

void
_acb_pow_arb_exp(acb_t z, const acb_t x, const arb_t y, long prec)
{
    acb_t t;
    acb_init(t);
    acb_log(t, x, prec);
    acb_mul_arb(t, t, y, prec);
    acb_exp(z, t, prec);
    acb_clear(t);
}

#define BINEXP_LIMIT 64

void
acb_pow_arb(acb_t z, const acb_t x, const arb_t y, long prec)
{
    const arf_struct * ymid = arb_midref(y);
    const mag_struct * yrad = arb_radref(y);

    if (arb_is_zero(y))
    {
        acb_one(z);
        return;
    }

    if (mag_is_zero(yrad))
    {
        /* small half-integer or integer */
        if (arf_cmpabs_2exp_si(ymid, BINEXP_LIMIT) < 0 &&
            arf_is_int_2exp_si(ymid, -1))
        {
            fmpz_t e;
            fmpz_init(e);            

            if (arf_is_int(ymid))
            {
                arf_get_fmpz_fixed_si(e, ymid, 0);
                acb_pow_fmpz_binexp(z, x, e, prec);
            }
            else
            {
                /* hack: give something finite here (should fix sqrt/rsqrt etc) */
                if (arb_contains_zero(acb_imagref(x)) && arb_contains_nonpositive(acb_realref(x)))
                {
                    _acb_pow_arb_exp(z, x, y, prec);
                    fmpz_clear(e);
                    return;
                }

                arf_get_fmpz_fixed_si(e, ymid, -1);
                acb_sqrt(z, x, prec + fmpz_bits(e));
                acb_pow_fmpz_binexp(z, z, e, prec);
            }

            fmpz_clear(e);
            return;
        }
    }

    _acb_pow_arb_exp(z, x, y, prec);
}

void
acb_pow(acb_t z, const acb_t x, const acb_t y, long prec)
{
    if (arb_is_zero(acb_imagref(y)))
    {
        acb_pow_arb(z, x, acb_realref(y), prec);
    }
    else
    {
        _acb_pow_exp(z, x, y, prec);
    }
}

