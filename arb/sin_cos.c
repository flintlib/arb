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

    Copyright (C) 2012-2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

#define MAGLIM(prec) FLINT_MAX(65536, (4*prec))

static void
_arf_sin(arf_t z, const arf_t x, long prec, arf_rnd_t rnd)
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

    mpfr_sin(zf, xf, arf_rnd_to_mpfr(rnd));

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

static void
_arf_cos(arf_t z, const arf_t x, long prec, arf_rnd_t rnd)
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

    mpfr_cos(zf, xf, arf_rnd_to_mpfr(rnd));

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

static void
_arf_sin_cos(arf_t z, arf_t w, const arf_t x, long prec, arf_rnd_t rnd)
{
    mpfr_t xf, zf, wf;
    mp_ptr zptr, wptr, tmp, tmp2;
    mp_srcptr xptr;
    mp_size_t xn, zn, wn, val;
    TMP_INIT;
    TMP_START;

    zn = wn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    tmp = TMP_ALLOC(2 * zn * sizeof(mp_limb_t));
    tmp2 = tmp + zn;

    ARF_GET_MPN_READONLY(xptr, xn, x);

    xf->_mpfr_d = (mp_ptr) xptr;
    xf->_mpfr_prec = xn * FLINT_BITS;
    xf->_mpfr_sign = ARF_SGNBIT(x) ? -1 : 1;
    xf->_mpfr_exp = ARF_EXP(x);

    zf->_mpfr_d = tmp;
    zf->_mpfr_prec = prec;
    zf->_mpfr_sign = 1;
    zf->_mpfr_exp = 0;

    wf->_mpfr_d = tmp2;
    wf->_mpfr_prec = prec;
    wf->_mpfr_sign = 1;
    wf->_mpfr_exp = 0;

    mpfr_sin_cos(zf, wf, xf, arf_rnd_to_mpfr(rnd));

    val = 0;
    while (tmp[val] == 0)
        val++;
    ARF_GET_MPN_WRITE(zptr, zn - val, z);
    flint_mpn_copyi(zptr, tmp + val, zn - val);
    if (zf->_mpfr_sign < 0)
        ARF_NEG(z);
    fmpz_set_si(ARF_EXPREF(z), zf->_mpfr_exp);

    val = 0;
    while (tmp2[val] == 0)
        val++;
    ARF_GET_MPN_WRITE(wptr, wn - val, w);
    flint_mpn_copyi(wptr, tmp2 + val, wn - val);
    if (wf->_mpfr_sign < 0)
        ARF_NEG(w);
    fmpz_set_si(ARF_EXPREF(w), wf->_mpfr_exp);

    TMP_END;
}

void
arb_sin_arf(arb_t s, const arf_t x, long prec, long maglim)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            arb_zero(s);
        }
        else if (arf_is_nan(x))
        {
            arb_indeterminate(s);
        }
        else
        {
            arf_zero(arb_midref(s));
            mag_one(arb_radref(s));
        }
    }
    else
    {
        long xmag;

        /* 2^(xmag-1) <= |x| < 2^xmag */
        xmag = ARF_EXP(x);

        if (xmag >= -(prec/3) - 2 && xmag <= maglim)
        {
            _arf_sin(arb_midref(s), x, prec, ARF_RND_DOWN);
            /* must be inexact */
            arf_mag_set_ulp(arb_radref(s), arb_midref(s), prec);
        }
        /* sin x = x + eps, |eps| < x^3 */
        else if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_mul_ui(t, ARF_EXPREF(x), 3);
            arb_set_arf(s, x);
            arb_set_round(s, s, prec);
            arb_add_error_2exp_fmpz(s, t);
            fmpz_clear(t);
        }
        /* huge */
        else
        {
            arf_zero(arb_midref(s));
            mag_one(arb_radref(s));
        }
    }
}

void
arb_cos_arf(arb_t c, const arf_t x, long prec, long maglim)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            arb_one(c);
        }
        else if (arf_is_nan(x))
        {
            arb_indeterminate(c);
        }
        else
        {
            arf_zero(arb_midref(c));
            mag_one(arb_radref(c));
        }
    }
    else
    {
        long xmag;

        /* 2^(xmag-1) <= |x| < 2^xmag */
        xmag = ARF_EXP(x);

        if (xmag >= -(prec/2) - 2 && xmag <= maglim)
        {
            _arf_cos(arb_midref(c), x, prec, ARF_RND_DOWN);
            /* must be inexact */
            arf_mag_set_ulp(arb_radref(c), arb_midref(c), prec);
        }
        /* cos x = 1 - eps, |eps| < x^2 */
        else if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_mul_ui(t, ARF_EXPREF(x), 2);
            arb_one(c);
            arb_add_error_2exp_fmpz(c, t);
            fmpz_clear(t);
        }
        /* huge */
        else
        {
            arf_zero(arb_midref(c));
            mag_one(arb_radref(c));
        }
    }
}

void
arb_sin_cos_arf(arb_t s, arb_t c, const arf_t x, long prec, long maglim)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            arb_zero(s);
            arb_one(c);
        }
        else if (arf_is_nan(x))
        {
            arb_indeterminate(s);
            arb_set(c, s);
        }
        else
        {
            arf_zero(arb_midref(s));
            mag_one(arb_radref(s));
            arb_set(c, s);
        }
    }
    else
    {
        long xmag;

        /* 2^(xmag-1) <= |x| < 2^xmag */
        xmag = ARF_EXP(x);

        if (xmag >= -(prec/2) - 2 && xmag <= maglim)
        {
            _arf_sin_cos(arb_midref(s), arb_midref(c), x, prec, ARF_RND_DOWN);
            /* must be inexact */
            arf_mag_set_ulp(arb_radref(s), arb_midref(s), prec);
            arf_mag_set_ulp(arb_radref(c), arb_midref(c), prec);
        }
        /* sin x = x + eps, |eps| < x^3 */
        /* cos x = 1 - eps, |eps| < x^2 */
        else if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_mul_ui(t, ARF_EXPREF(x), 3);
            arb_set_arf(s, x);
            arb_set_round(s, s, prec);
            arb_add_error_2exp_fmpz(s, t);
            fmpz_divexact_ui(t, t, 3);
            fmpz_mul_ui(t, t, 2);
            arb_one(c);
            arb_add_error_2exp_fmpz(c, t);
            fmpz_clear(t);
        }
        /* huge */
        else
        {
            arf_zero(arb_midref(s));
            mag_one(arb_radref(s));
            arb_set(c, s);
        }
    }
}

void
arb_sin(arb_t s, const arb_t x, long prec)
{
    if (arb_is_exact(x))
    {
        arb_sin_arf(s, arb_midref(x), prec, MAGLIM(prec));
    }
    else
    {
        mag_t t;
        mag_init(t);

        if (mag_cmp_2exp_si(arb_radref(x), 1) > 0)
            mag_set_ui_2exp_si(t, 1, 1);
        else
            mag_set(t, arb_radref(x));

        arb_sin_arf(s, arb_midref(x), prec, MAGLIM(prec));
        mag_add(arb_radref(s), arb_radref(s), t);

        mag_clear(t);
    }
}

void
arb_cos(arb_t s, const arb_t x, long prec)
{
    if (arb_is_exact(x))
    {
        arb_cos_arf(s, arb_midref(x), prec, MAGLIM(prec));
    }
    else
    {
        mag_t t;
        mag_init(t);

        if (mag_cmp_2exp_si(arb_radref(x), 1) > 0)
            mag_set_ui_2exp_si(t, 1, 1);
        else
            mag_set(t, arb_radref(x));

        arb_cos_arf(s, arb_midref(x), prec, MAGLIM(prec));
        mag_add(arb_radref(s), arb_radref(s), t);

        mag_clear(t);
    }
}

void
arb_sin_cos(arb_t s, arb_t c, const arb_t x, long prec)
{
    if (arb_is_exact(x))
    {
        arb_sin_cos_arf(s, c, arb_midref(x), prec, MAGLIM(prec));
    }
    else
    {
        mag_t t;
        mag_init(t);

        if (mag_cmp_2exp_si(arb_radref(x), 1) > 0)
            mag_set_ui_2exp_si(t, 1, 1);
        else
            mag_set(t, arb_radref(x));

        arb_sin_cos_arf(s, c, arb_midref(x), prec, MAGLIM(prec));
        mag_add(arb_radref(s), arb_radref(s), t);
        mag_add(arb_radref(c), arb_radref(c), t);

        mag_clear(t);
    }
}

