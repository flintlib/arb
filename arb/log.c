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

#include "arb.h"

#define BIG_EXPONENT_BITS 20
#define BIG_EXPONENT (1L << BIG_EXPONENT_BITS)

/* requires x != 1 */
static void
_arf_log(arf_t z, const arf_t x, long prec, arf_rnd_t rnd)
{
    mpfr_t xf, zf;
    mp_ptr zptr, tmp;
    mp_srcptr xptr;
    mp_size_t xn, zn, val;
    TMP_INIT;
    TMP_START;

    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    tmp = TMP_ALLOC(zn * sizeof(mp_limb_t));

    ARF_GET_MPN_READONLY(xptr, xn, x);

    xf->_mpfr_d = (mp_ptr) xptr;
    xf->_mpfr_prec = xn * FLINT_BITS;
    xf->_mpfr_sign = ARF_SGNBIT(x) ? -1 : 1;
    xf->_mpfr_exp = ARF_EXP(x);

    zf->_mpfr_d = tmp;
    zf->_mpfr_prec = prec;
    zf->_mpfr_sign = 1;
    zf->_mpfr_exp = 0;

    mpfr_log(zf, xf, arf_rnd_to_mpfr(rnd));

    val = 0;
    while (tmp[val] == 0)
        val++;

    ARF_GET_MPN_WRITE(zptr, zn - val, z);
    flint_mpn_copyi(zptr, tmp + val, zn - val);
    if (zf->_mpfr_sign < 0)
        ARF_NEG(z);

    fmpz_set_si(ARF_EXPREF(z), zf->_mpfr_exp);

    TMP_END;
}

void
arb_log_arf(arb_t y, const arf_t x, long prec)
{
    fmpz_t exp;

    if (arf_is_special(x))
    {
        if (arf_is_pos_inf(x))
            arb_pos_inf(y);
        else
            arb_indeterminate(y);
        return;
    }

    if (ARF_SGNBIT(x))
    {
        arb_indeterminate(y);
        return;
    }

    if (ARF_IS_POW2(x))
    {
        if (fmpz_is_one(ARF_EXPREF(x)))
        {
            arb_zero(y);
        }
        else
        {
            fmpz_init(exp);
            _fmpz_add_fast(exp, ARF_EXPREF(x), -1);
            arb_const_log2(y, prec + 2);
            arb_mul_fmpz(y, y, exp, prec);
            fmpz_clear(exp);
        }
        return;
    }

    if (ARF_EXP(x) >= -BIG_EXPONENT && ARF_EXP(x) <= BIG_EXPONENT)
    {
        _arf_log(arb_midref(y), x, prec, ARB_RND);
        arf_mag_set_ulp(arb_radref(y), arb_midref(y), prec);
    }
    else
    {
        arf_t t;
        arb_t c;
        long wp;

        arf_init(t);
        arb_init(c);
        fmpz_init(exp);

        fmpz_neg(exp, ARF_EXPREF(x));
        arf_mul_2exp_fmpz(t, x, exp);

        wp = prec + 4 - fmpz_bits(exp);
        wp = FLINT_MAX(wp, 4);

        _arf_log(arb_midref(y), t, wp, ARB_RND);
        arf_mag_set_ulp(arb_radref(y), arb_midref(y), wp);

        arb_const_log2(c, prec + 4);
        arb_mul_fmpz(c, c, exp, prec + 4);

        arb_sub(y, y, c, prec);

        arf_clear(t);
        arb_clear(c);
        fmpz_clear(exp);
    }
}

void
arb_log(arb_t y, const arb_t x, long prec)
{
    if (arb_is_exact(x))
    {
        arb_log_arf(y, arb_midref(x), prec);
    }
    else
    {
        /*
        Let the input be [a-b, a+b]. We require a > b >= 0 (otherwise the
        interval contains zero or a negative number and the logarithm is not
        defined). The error is largest at a-b, and we have

        log(a) - log(a-b) = log(1 + b/(a-b)).
        */
        mag_t err;
        mag_init(err);

        arb_get_mag_infimum_lower(err, x);

        if (mag_is_zero(err))
        {
            mag_inf(err);
        }
        else
        {
            mag_div(err, arb_radref(x), err);
            mag_log1p(err, err);
        }

        arb_log_arf(y, arb_midref(x), prec);

        mag_add(arb_radref(y), arb_radref(y), err);
        mag_clear(err);
    }
}

void
arb_log_ui(arb_t z, ulong x, long prec)
{
    arf_t t;
    arf_init(t);
    arf_set_ui(t, x);
    arb_log_arf(z, t, prec);
    arf_clear(t);
}

void
arb_log_fmpz(arb_t z, const fmpz_t x, long prec)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, x);
    arb_log_arf(z, t, prec);
    arf_clear(t);
}
