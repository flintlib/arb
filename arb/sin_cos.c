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

    Copyright (C) 2012-2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"
#include "mpn_extras.h"

#define TMP_ALLOC_LIMBS(__n) TMP_ALLOC((__n) * sizeof(mp_limb_t))

int _arf_get_integer_mpn(mp_ptr y, mp_srcptr x, mp_size_t xn, long exp);

int _arf_set_mpn_fixed(arf_t z, mp_srcptr xp, mp_size_t xn, mp_size_t fixn, int negative, long prec);

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

    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);

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

    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);

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

    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);

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
arb_sin_cos_arf_new(arb_t zsin, arb_t zcos, const arf_t x, long prec)
{
    int want_sin, want_cos;
    long exp, wp, wn, N, r, wprounded;
    mp_ptr tmp, w, sina, cosa, sinb, cosb, ta, tb;
    mp_ptr sinptr, cosptr;
    mp_limb_t p1, q1bits, p2, q2bits, error, error2;
    int negative, inexact, octant;
    int sinnegative, cosnegative, swapsincos;
    TMP_INIT;

    want_sin = (zsin != NULL);
    want_cos = (zcos != NULL);

    exp = ARF_EXP(x);
    negative = ARF_SGNBIT(x);

    /* Absolute working precision (NOT rounded to a limb multiple) */
    wp = prec + 8;
    if (want_sin && exp <= 0)
        wp += (-exp);
    /* Number of limbs */
    wn = (wp + FLINT_BITS - 1) / FLINT_BITS;
    /* Precision rounded to a number of bits */
    wprounded = FLINT_BITS * wn;
    /* Don't be close to the boundary (to allow adding adding the
       Taylor series truncation error without overflow) */
    wp = FLINT_MAX(wp, wprounded - (FLINT_BITS - 4));

    /* Too high precision to use table -- use generic algorithm */
    if (wp > ARB_SIN_COS_TAB2_PREC)
    {
        if (want_sin && want_cos)
        {
            _arf_sin_cos(arb_midref(zsin), arb_midref(zcos), x, prec, ARB_RND);
            arf_mag_set_ulp(arb_radref(zsin), arb_midref(zsin), prec);
            arf_mag_set_ulp(arb_radref(zcos), arb_midref(zcos), prec);
        }
        else if (want_sin)
        {
            _arf_sin(arb_midref(zsin), x, prec, ARB_RND);
            arf_mag_set_ulp(arb_radref(zsin), arb_midref(zsin), prec);
        }
        else
        {
            _arf_cos(arb_midref(zcos), x, prec, ARB_RND);
            arf_mag_set_ulp(arb_radref(zcos), arb_midref(zcos), prec);
        }
        return;
    }

    TMP_START;
    tmp = TMP_ALLOC_LIMBS(9 * wn);
    w    = tmp;         /* requires wn limbs */
    sina = w    + wn;   /* requires wn limbs */
    cosa = sina + wn;   /* requires wn limbs */
    sinb = cosa + wn;   /* requires wn limbs */
    cosb = sinb + wn;   /* requires wn limbs */
    ta   = cosb + wn;   /* requires 2*wn limbs */
    tb   = ta + 2*wn;   /* requires 2*wn limbs */

    /* reduce modulo pi/4 */
    if (_arb_get_mpn_fixed_mod_pi4(w, NULL, &octant, &error, x, wn) == 0)
    {
        /* may run out of precision for pi/4 */
        if (want_sin && want_cos)
        {
            _arf_sin_cos(arb_midref(zsin), arb_midref(zcos), x, prec, ARB_RND);
            arf_mag_set_ulp(arb_radref(zsin), arb_midref(zsin), prec);
            arf_mag_set_ulp(arb_radref(zcos), arb_midref(zcos), prec);
        }
        else if (want_sin)
        {
            _arf_sin(arb_midref(zsin), x, prec, ARB_RND);
            arf_mag_set_ulp(arb_radref(zsin), arb_midref(zsin), prec);
        }
        else
        {
            _arf_cos(arb_midref(zcos), x, prec, ARB_RND);
            arf_mag_set_ulp(arb_radref(zcos), arb_midref(zcos), prec);
        }
        TMP_END;
        return;
    }

    sinnegative = (octant >= 4) ^ negative;
    cosnegative = (octant >= 2 && octant <= 5);
    swapsincos = (octant == 1 || octant == 2 || octant == 5 || octant == 6);

    /* Table-based argument reduction (1 or 2 steps) */
    if (wp <= ARB_SIN_COS_TAB1_PREC)
    {
        q1bits = ARB_SIN_COS_TAB1_BITS;
        q2bits = 0;

        p1 = w[wn-1] >> (FLINT_BITS - q1bits);
        w[wn-1] -= (p1 << (FLINT_BITS - q1bits));
        p2 = 0;
    }
    else
    {
        q1bits = ARB_SIN_COS_TAB21_BITS;
        q2bits = ARB_SIN_COS_TAB21_BITS + ARB_SIN_COS_TAB22_BITS;

        p1 = w[wn-1] >> (FLINT_BITS - q1bits);
        w[wn-1] -= (p1 << (FLINT_BITS - q1bits));
        p2 = w[wn-1] >> (FLINT_BITS - q2bits);
        w[wn-1] -= (p2 << (FLINT_BITS - q2bits));
    }

    /* |w| <= 2^-r */
    r = _arb_mpn_leading_zeros(w, wn);

    /* Choose number of terms N such that Taylor series truncation
       error is <= 2^-wp */
    N = _arb_exp_taylor_bound(-r, wp);

    /* the summation for sin/cos is actually done to (2N-1)! */
    N = (N + 1) / 2;

    if (N < 14)
    {
        /* Evaluate Taylor series */
        _arb_sin_cos_taylor_rs(sina, cosa, &error2, w, wn, N, 0, 1);
        /* Taylor series evaluation error */
        error += error2;
        /* Taylor series truncation error */
        error += 1UL << (wprounded-wp);
    }
    else  /* Compute cos(a) from sin(a) using a square root. */
    {
        /* Evaluate Taylor series */
        _arb_sin_cos_taylor_rs(sina, cosa, &error2, w, wn, N, 1, 1);
        error += error2;
        error += 1UL << (wprounded-wp);

        if (flint_mpn_zero_p(sina, wn))
        {
            flint_mpn_store(cosa, wn, LIMB_ONES);
            error = FLINT_MAX(error, 1);
        }
        else
        {
            mpn_sqr(ta, sina, wn);
            /* 1 - s^2 (negation guaranteed to have borrow) */
            mpn_neg(ta, ta, 2 * wn);
            /* top limb of ta must be nonzero, so no need to normalize
               before calling mpn_sqrtrem */
            mpn_sqrtrem(cosa, ta, ta, 2 * wn);

            /* When converting sin to cos, the error for cos must be
               smaller than the error for sin; but we also get 1 ulp
               extra error from the square root. */
            error += 1;
        }
    }

    /*
    sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
    cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)

    [F1+e1]*[F2+e2] + [F3+e3]*[F4+e4] - F1 F2 + F3 F4 = 

       e1 e2 + e3 e4 + e2 F1 + e1 F2 + e4 F3 + e3 F4

       <= (e1 + e2 + e3 + e4) + 1    (ulp)

       <= 2*left_err + 2*right_err + 1     (ulp)

    Truncating both terms before adding adds another 2 ulp, so the error
    is bounded by 

       <= 2*left_err + 2*right_err + 3     (ulp)
    */

    if (p1 == 0 && p2 == 0)  /* no table lookups */
    {
        sinptr = sina;
        cosptr = cosa;
    }
    else if (p1 == 0 || p2 == 0)    /* only one table lookup */
    {
        mp_srcptr sinc, cosc;

        if (wp <= ARB_SIN_COS_TAB1_PREC)  /* must be in table 1 */
        {
            sinc = arb_sin_cos_tab1[2 * p1] + ARB_SIN_COS_TAB1_LIMBS - wn;
            cosc = arb_sin_cos_tab1[2 * p1 + 1] + ARB_SIN_COS_TAB1_LIMBS - wn;
        }
        else if (p1 != 0)
        {
            sinc = arb_sin_cos_tab21[2 * p1] + ARB_SIN_COS_TAB2_LIMBS - wn;
            cosc = arb_sin_cos_tab21[2 * p1 + 1] + ARB_SIN_COS_TAB2_LIMBS - wn;
        }
        else
        {
            sinc = arb_sin_cos_tab22[2 * p2] + ARB_SIN_COS_TAB2_LIMBS - wn;
            cosc = arb_sin_cos_tab22[2 * p2 + 1] + ARB_SIN_COS_TAB2_LIMBS - wn;
        }

        if ((want_sin && !swapsincos) || (want_cos && swapsincos))
        {
            mpn_mul_n(ta, sina, cosc, wn);
            mpn_mul_n(tb, cosa, sinc, wn);
            mpn_add_n(w, ta + wn, tb + wn, wn);
        }

        if ((want_cos && !swapsincos) || (want_sin && swapsincos))
        {
            mpn_mul_n(ta, cosa, cosc, wn);
            mpn_mul_n(tb, sina, sinc, wn);
            mpn_sub_n(ta, ta + wn, tb + wn, wn);
        }

        sinptr = w;
        cosptr = ta;

        error = 2 * error + 2 * 1 + 3;
    }
    else        /* two table lookups, must be in table 2 */
    {
        mp_srcptr sinc, cosc, sind, cosd;

        sinc = arb_sin_cos_tab21[2 * p1] + ARB_SIN_COS_TAB2_LIMBS - wn;
        cosc = arb_sin_cos_tab21[2 * p1 + 1] + ARB_SIN_COS_TAB2_LIMBS - wn;
        sind = arb_sin_cos_tab22[2 * p2] + ARB_SIN_COS_TAB2_LIMBS - wn;
        cosd = arb_sin_cos_tab22[2 * p2 + 1] + ARB_SIN_COS_TAB2_LIMBS - wn;

        mpn_mul_n(ta, sinc, cosd, wn);
        mpn_mul_n(tb, cosc, sind, wn);
        mpn_add_n(sinb, ta + wn, tb + wn, wn);

        mpn_mul_n(ta, cosc, cosd, wn);
        mpn_mul_n(tb, sinc, sind, wn);
        mpn_sub_n(cosb, ta + wn, tb + wn, wn);

        error2 = 2 * 1 + 2 * 1 + 3;

        if ((want_sin && !swapsincos) || (want_cos && swapsincos))
        {
            mpn_mul_n(ta, sina, cosb, wn);
            mpn_mul_n(tb, cosa, sinb, wn);
            mpn_add_n(w, ta + wn, tb + wn, wn);
        }

        if ((want_cos && !swapsincos) || (want_sin && swapsincos))
        {
            mpn_mul_n(ta, cosa, cosb, wn);
            mpn_mul_n(tb, sina, sinb, wn);
            mpn_sub_n(ta, ta + wn, tb + wn, wn);
        }

        error = 2 * error + 2 * error2 + 3;

        sinptr = w;
        cosptr = ta;
    }

    if (swapsincos)
    {
        mp_ptr tmptr = sinptr;
        sinptr = cosptr;
        cosptr = tmptr;
    }

    /* The accumulated error */
    if (want_sin)
    {
        mag_set_ui_2exp_si(arb_radref(zsin), error, -wprounded);

        if (want_cos)
            mag_set(arb_radref(zcos), arb_radref(zsin));
    }
    else
    {
        mag_set_ui_2exp_si(arb_radref(zcos), error, -wprounded);
    }

    /* Set the midpoint */
    if (want_sin)
    {
        inexact = _arf_set_mpn_fixed(arb_midref(zsin), sinptr,
            wn, wn, sinnegative, prec);
        if (inexact)
            arf_mag_add_ulp(arb_radref(zsin),
                arb_radref(zsin), arb_midref(zsin), prec);
    }

    if (want_cos)
    {
        inexact = _arf_set_mpn_fixed(arb_midref(zcos), cosptr,
            wn, wn, cosnegative, prec);
        if (inexact)
            arf_mag_add_ulp(arb_radref(zcos),
                arb_radref(zcos), arb_midref(zcos), prec);
    }

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
            arb_sin_cos_arf_new(s, NULL, x, prec);
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
            arb_sin_cos_arf_new(NULL, c, x, prec);
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
            arb_sin_cos_arf_new(s, c, x, prec);
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

