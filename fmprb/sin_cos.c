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

#include "fmprb.h"

#define MAGLIM(prec) FLINT_MAX(65536, (4*prec))

static long
_fmpr_sin(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    long r;
    CALL_MPFR_FUNC(r, mpfr_sin, y, x, prec, rnd);
    return r;
}

static long
_fmpr_cos(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    long r;
    CALL_MPFR_FUNC(r, mpfr_cos, y, x, prec, rnd);
    return r;
}

static void
_fmpr_sin_cos(long * r1, long * r2, fmpr_t s, fmpr_t c, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    CALL_MPFR_FUNC_2X1(*r1, *r2, mpfr_sin_cos, s, c, x, prec, rnd);
}

void
fmprb_sin_fmpr(fmprb_t s, const fmpr_t x, long prec, long maglim)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
        {
            fmprb_zero(s);
        }
        else if (fmpr_is_nan(x))
        {
            fmpr_nan(fmprb_midref(s));
            fmpr_pos_inf(fmprb_radref(s));
        }
        else
        {
            fmpr_zero(fmprb_midref(s));
            fmpr_one(fmprb_radref(s));
        }
    }
    else
    {
        long r;
        fmpz_t mag;

        fmpz_init(mag);

        /* 2^(mag-1) <= |x| < 2^mag */
        fmpz_add_ui(mag, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));

        /* sin x = x + eps, |eps| < x^3 */
        if (fmpz_cmp_si(mag, -(prec/3) - 2) < 0)
        {
            fmpz_mul_ui(mag, mag, 3);
            fmprb_set_fmpr(s, x);
            fmprb_set_round(s, s, prec);
            fmprb_add_error_2exp_fmpz(s, mag);
        }
        /* huge */
        else if (fmpz_cmp_si(mag, maglim) > 0)
        {
            fmpr_zero(fmprb_midref(s));
            fmpr_one(fmprb_radref(s));
        }
        else
        {
            r = _fmpr_sin(fmprb_midref(s), x, prec, FMPR_RND_DOWN);
            fmpr_set_error_result(fmprb_radref(s), fmprb_midref(s), r);
        }

        fmpz_clear(mag);
    }
}

void
fmprb_cos_fmpr(fmprb_t c, const fmpr_t x, long prec, long maglim)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
        {
            fmprb_one(c);
        }
        else if (fmpr_is_nan(x))
        {
            fmpr_nan(fmprb_midref(c));
            fmpr_pos_inf(fmprb_radref(c));
        }
        else
        {
            fmpr_zero(fmprb_midref(c));
            fmpr_one(fmprb_radref(c));
        }
    }
    else
    {
        long r;
        fmpz_t mag;

        fmpz_init(mag);

        /* 2^(mag-1) <= |x| < 2^mag */
        fmpz_add_ui(mag, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));

        /* cos x = 1 - eps, |eps| < x^2 */
        if (fmpz_cmp_si(mag, -(prec/2) - 2) < 0)
        {
            fmpz_mul_ui(mag, mag, 2);
            fmprb_one(c);
            fmprb_add_error_2exp_fmpz(c, mag);
        }
        /* huge */
        else if (fmpz_cmp_si(mag, maglim) > 0)
        {
            fmpr_zero(fmprb_midref(c));
            fmpr_one(fmprb_radref(c));
        }
        else
        {
            r = _fmpr_cos(fmprb_midref(c), x, prec, FMPR_RND_DOWN);
            fmpr_set_error_result(fmprb_radref(c), fmprb_midref(c), r);
        }

        fmpz_clear(mag);
    }
}

void
fmprb_sin_cos_fmpr(fmprb_t s, fmprb_t c, const fmpr_t x, long prec, long maglim)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
        {
            fmprb_zero(s);
            fmprb_one(c);
        }
        else if (fmpr_is_nan(x))
        {
            fmpr_nan(fmprb_midref(s));
            fmpr_pos_inf(fmprb_radref(s));
            fmprb_set(c, s);
        }
        else
        {
            fmpr_zero(fmprb_midref(s));
            fmpr_one(fmprb_radref(s));
            fmprb_set(c, s);
        }
    }
    else
    {
        long r1, r2;
        fmpz_t mag;

        fmpz_init(mag);

        /* 2^(mag-1) <= |x| < 2^mag */
        fmpz_add_ui(mag, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));

        /* sin x = x + eps, |eps| < x^3 */
        /* cos x = 1 - eps, |eps| < x^2 */
        if (fmpz_cmp_si(mag, -(prec/2) - 2) < 0)
        {
            fmpz_mul_ui(mag, mag, 3);
            fmprb_set_fmpr(s, x);
            fmprb_set_round(s, s, prec);
            fmprb_add_error_2exp_fmpz(s, mag);
            fmpz_divexact_ui(mag, mag, 3);
            fmpz_mul_ui(mag, mag, 2);
            fmprb_one(c);
            fmprb_add_error_2exp_fmpz(c, mag);
        }
        /* huge */
        else if (fmpz_cmp_si(mag, maglim) > 0)
        {
            fmpr_zero(fmprb_midref(s));
            fmpr_one(fmprb_radref(s));
            fmprb_set(c, s);
        }
        else
        {
            _fmpr_sin_cos(&r1, &r2, fmprb_midref(s), fmprb_midref(c), x, prec, FMPR_RND_DOWN);
            fmpr_set_error_result(fmprb_radref(s), fmprb_midref(s), r1);
            fmpr_set_error_result(fmprb_radref(c), fmprb_midref(c), r2);
        }

        fmpz_clear(mag);
    }
}

void
fmprb_sin(fmprb_t s, const fmprb_t x, long prec)
{
    if (fmprb_is_exact(x))
    {
        fmprb_sin_fmpr(s, fmprb_midref(x), prec, MAGLIM(prec));
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);
        if (fmpr_cmpabs_2exp_si(fmprb_radref(x), 1) > 0)
            fmpr_set_si_2exp_si(t, 1, 1);
        else
            fmpr_set(t, fmprb_radref(x));
        fmprb_sin_fmpr(s, fmprb_midref(x), prec, MAGLIM(prec));
        fmpr_add(fmprb_radref(s), fmprb_radref(s), t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_clear(t);
    }
    fmprb_adjust(s);
}

void
fmprb_cos(fmprb_t c, const fmprb_t x, long prec)
{
    if (fmprb_is_exact(x))
    {
        fmprb_cos_fmpr(c, fmprb_midref(x), prec, MAGLIM(prec));
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);
        if (fmpr_cmpabs_2exp_si(fmprb_radref(x), 1) > 0)
            fmpr_set_si_2exp_si(t, 1, 1);
        else
            fmpr_set(t, fmprb_radref(x));
        fmprb_cos_fmpr(c, fmprb_midref(x), prec, MAGLIM(prec));
        fmpr_add(fmprb_radref(c), fmprb_radref(c), t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_clear(t);
    }
    fmprb_adjust(c);
}

void
fmprb_sin_cos(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)
{
    if (fmprb_is_exact(x))
    {
        fmprb_sin_cos_fmpr(s, c, fmprb_midref(x), prec, MAGLIM(prec));
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);
        if (fmpr_cmpabs_2exp_si(fmprb_radref(x), 1) > 0)
            fmpr_set_si_2exp_si(t, 1, 1);
        else
            fmpr_set(t, fmprb_radref(x));
        fmprb_sin_cos_fmpr(s, c, fmprb_midref(x), prec, MAGLIM(prec));
        fmpr_add(fmprb_radref(s), fmprb_radref(s), t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_add(fmprb_radref(c), fmprb_radref(c), t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_clear(t);
    }
    fmprb_adjust(s);
    fmprb_adjust(c);
}

