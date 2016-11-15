/*
    Copyright (C) 2012-2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/mpn_extras.h"
#include "arb.h"

#define TMP_ALLOC_LIMBS(__n) TMP_ALLOC((__n) * sizeof(mp_limb_t))
#define MAGLIM(prec) FLINT_MAX(65536, 4*prec)

static void
_arf_sin(arf_t z, const arf_t x, slong prec, arf_rnd_t rnd)
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
_arf_cos(arf_t z, const arf_t x, slong prec, arf_rnd_t rnd)
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
_arf_sin_cos(arf_t z, arf_t w, const arf_t x, slong prec, arf_rnd_t rnd)
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
arb_zero_pm_one(arb_t res)
{
    arf_zero(arb_midref(res));
    mag_one(arb_radref(res));
}

static __inline__ void
mag_nonzero_fast_mul(mag_t z, const mag_t x, const mag_t y)
{
    MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) + LIMB_ONE;
    MAG_EXP(z) = MAG_EXP(x) + MAG_EXP(y);
    MAG_FAST_ADJUST_ONE_TOO_SMALL(z);
}

static __inline__ void
mag_nonzero_fast_add(mag_t z, const mag_t x, const mag_t y)
{
    slong shift = MAG_EXP(x) - MAG_EXP(y);

    if (shift == 0)
    {
        MAG_EXP(z) = MAG_EXP(x);
        MAG_MAN(z) = MAG_MAN(x) + MAG_MAN(y);
        MAG_FAST_ADJUST_ONE_TOO_LARGE(z); /* may need two adjustments */
    }
    else if (shift > 0)
    {
        MAG_EXP(z) = MAG_EXP(x);

        if (shift >= MAG_BITS)
            MAG_MAN(z) = MAG_MAN(x) + LIMB_ONE;
        else
            MAG_MAN(z) = MAG_MAN(x) + (MAG_MAN(y) >> shift) + LIMB_ONE;
    }
    else
    {
        shift = -shift;
        MAG_EXP(z) = MAG_EXP(y);

        if (shift >= MAG_BITS)
            MAG_MAN(z) = MAG_MAN(y) + LIMB_ONE;
        else
            MAG_MAN(z) = MAG_MAN(y) + (MAG_MAN(x) >> shift) + LIMB_ONE;
    }

    MAG_FAST_ADJUST_ONE_TOO_LARGE(z);
}

static __inline__ int
mag_nonzero_fast_cmp(const mag_t x, const mag_t y)
{
    if (MAG_EXP(x) == MAG_EXP(y))
        return (MAG_MAN(x) < MAG_MAN(y)) ? -1 : 1;
    else
        return (MAG_EXP(x) < MAG_EXP(y)) ? -1 : 1;
}

static __inline__ void
mag_fast_set(mag_t x, const mag_t y)
{
    MAG_EXP(x) = MAG_EXP(y);
    MAG_MAN(x) = MAG_MAN(y);
}

void
arb_sin_cos_tiny_huge(arb_t s, arb_t c, const arf_t x, const mag_t xrad, slong prec)
{
    /* sin(x) = x + eps, |eps| <= x^3/6 */
    /* cos(x) = 1 - eps, |eps| <= x^2/2 */
    if (fmpz_sgn(ARF_EXPREF(x)) < 0)
    {
        mag_t t, u;
        mag_init(t);
        mag_init(u);

        /* TODO: min by 2 */

        arf_get_mag(t, x);
        mag_add(t, t, xrad);
        mag_mul(u, t, t);

        if (s != NULL)
        {
            mag_mul(t, u, t);
            mag_div_ui(t, t, 6);
            mag_add(t, t, xrad);

            arb_set_arf(s, x);
            arb_set_round(s, s, prec);
            arb_add_error_mag(s, t);
        }

        if (c != NULL)
        {
            arb_one(c);
            mag_mul_2exp_si(arb_radref(c), u, -1);
        }

        mag_clear(t);
        mag_clear(u);
    }
    else
    {
        if (s != NULL) arb_zero_pm_one(s);
        if (c != NULL) arb_zero_pm_one(c);
    }
}

void
arb_sin_cos_using_mpfr(arb_t s, arb_t c, const arf_t x, const mag_t xrad, slong prec)
{
    mag_t t;
    int want_sin, want_cos;

    want_sin = (s != NULL);
    want_cos = (c != NULL);

    mag_init(t);

    if (mag_cmp_2exp_si(xrad, 1) > 0)
        mag_set_ui_2exp_si(t, 1, 1);
    else
        mag_set(t, xrad);

    if (want_sin && want_cos)
    {
        _arf_sin_cos(arb_midref(s), arb_midref(c), x, prec, ARB_RND);
        arf_mag_set_ulp(arb_radref(s), arb_midref(s), prec);
        arf_mag_set_ulp(arb_radref(c), arb_midref(c), prec);
        mag_add(arb_radref(s), arb_radref(s), t);
        mag_add(arb_radref(c), arb_radref(c), t);
    }
    else if (want_sin)
    {
        _arf_sin(arb_midref(s), x, prec, ARB_RND);
        arf_mag_set_ulp(arb_radref(s), arb_midref(s), prec);
        mag_add(arb_radref(s), arb_radref(s), t);
    }
    else
    {
        _arf_cos(arb_midref(c), x, prec, ARB_RND);
        arf_mag_set_ulp(arb_radref(c), arb_midref(c), prec);
        mag_add(arb_radref(c), arb_radref(c), t);
    }

    mag_clear(t);
}

void
arb_sin_cos_special(arb_t s, arb_t c, const arf_t x, const mag_t xrad)
{
    int want_sin, want_cos;

    want_sin = (s != NULL);
    want_cos = (c != NULL);

    if (mag_is_inf(xrad))
    {
        if (arf_is_nan(x))
        {
            if (want_sin) arb_indeterminate(s);
            if (want_cos) arb_indeterminate(c);
        }
        else
        {
            if (want_sin) arb_zero_pm_one(s);
            if (want_cos) arb_zero_pm_one(c);
        }
        return;
    }

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            if (mag_is_zero(xrad))
            {
                if (want_sin) arb_zero(s);
                if (want_cos) arb_one(c);
            }
            else
            {
                mag_t t;
                mag_init_set(t, xrad);

                /* todo: could compute sin(+/- r), cos(+/- r) more accurately */
                if (mag_cmp_2exp_si(t, 1) < 0)
                {
                    if (want_sin)
                    {
                        arf_zero(arb_midref(s));
                        mag_set(arb_radref(s), t);
                    }

                    if (want_cos)
                    {
                        arf_one(arb_midref(c));
                        mag_mul(arb_radref(c), t, t);
                        mag_mul_2exp_si(arb_radref(c), arb_radref(c), -1);
                    }
                }
                else
                {
                    if (want_sin) arb_zero_pm_one(s);
                    if (want_cos) arb_zero_pm_one(c);
                }

                mag_clear(t);
            }
        }
        else
        {
            if (want_sin) arb_indeterminate(s);
            if (want_cos) arb_indeterminate(c);
        }
    }
}

void
arb_sin_cos_fast(arb_t zsin, arb_t zcos, const arf_t x, const mag_t xrad, slong prec)
{
    int want_sin, want_cos;
    slong radexp, exp, wp, wn, N, r, wprounded, maglim;
    mp_ptr tmp, w, sina, cosa, sinb, cosb, ta, tb;
    mp_ptr sinptr, cosptr;
    mp_limb_t p1, q1bits, p2, q2bits, error, error2, p1_tab1;
    int negative, inexact, octant;
    int sinnegative, cosnegative, swapsincos;
    mag_t xrad_copy;
    TMP_INIT;

    if (mag_is_inf(xrad) || arf_is_special(x))
    {
        arb_sin_cos_special(zsin, zcos, x, xrad);
        return;
    }

    want_sin = (zsin != NULL);
    want_cos = (zcos != NULL);

    /* From now on, both x and xrad are finite, and x is nonzero. */
    exp = ARF_EXP(x);
    negative = ARF_SGNBIT(x);
    maglim = MAGLIM(prec);
    radexp = MAG_EXP(xrad);

    /* Unlikely: tiny or huge midpoint (including any bignums). */
    if (exp < -(prec/2) - 2 || exp > maglim)
    {
        arb_sin_cos_tiny_huge(zsin, zcos, x, xrad, prec);
        return;
    }

    /* We will need a copy later in case of aliasing */
    *xrad_copy = *xrad;

    if (!mag_is_zero(xrad))
    {
        if (radexp <= 2 && radexp >= MAG_MIN_LAGOM_EXP)
        {
            /* Regular case: decrease precision to match generic max. accuracy. */
            /* Note: near x=0, the error is quadratic for cos. */
            if (want_cos && exp < -2)
                prec = FLINT_MIN(prec, 20 - 2 * radexp);
            else
                prec = FLINT_MIN(prec, 20 - radexp);
        }
        else if (fmpz_sgn(MAG_EXPREF(xrad)) > 0)  /* Huge radius. */
        {
            if (want_sin) arb_zero_pm_one(zsin);
            if (want_cos) arb_zero_pm_one(zcos);
            return;
        }
        else
        {
            /* The exponent could be too negative to use the fast
               bounds code, so set to an upper bound. */
            MAG_MAN(xrad_copy) = MAG_ONE_HALF;
            MAG_EXP(xrad_copy) = MAG_MIN_LAGOM_EXP + 1;
            /* important */
            radexp = MAG_EXP(xrad_copy);
        }
    }

    /* PART 2: the actual computation. */

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
        arb_sin_cos_using_mpfr(zsin, zcos, x, xrad, prec);
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
        arb_sin_cos_using_mpfr(zsin, zcos, x, xrad, prec);
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

        p1 = p1_tab1 = w[wn-1] >> (FLINT_BITS - q1bits);
        w[wn-1] -= (p1 << (FLINT_BITS - q1bits));
        p2 = 0;

        /* p1_tab1 will be used for the error bounds at the end. */
        p1_tab1 = p1;
    }
    else
    {
        q1bits = ARB_SIN_COS_TAB21_BITS;
        q2bits = ARB_SIN_COS_TAB21_BITS + ARB_SIN_COS_TAB22_BITS;

        /* p1_tab1 will be used for the error bounds at the end. */
        p1_tab1 = w[wn-1] >> (FLINT_BITS - ARB_SIN_COS_TAB1_BITS);

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
        error += UWORD(1) << (wprounded-wp);
    }
    else  /* Compute cos(a) from sin(a) using a square root. */
    {
        /* Evaluate Taylor series */
        _arb_sin_cos_taylor_rs(sina, cosa, &error2, w, wn, N, 1, 1);
        error += error2;
        error += UWORD(1) << (wprounded-wp);

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

    /* PART 3: compute propagated error and write output */

    if (swapsincos)
    {
        mp_ptr tmptr = sinptr;
        sinptr = cosptr;
        cosptr = tmptr;
    }

    /*
    We have two sources of error.

    1. Computation error:
        error * 2^(-wprounded)

    2. With input radius r != 0, the propagated error bound:
        sin(x):  min(2, r, |cos(x)|*r + 0.5*r^2)
        cos(x):  min(2, r, |sin(x)|*r + 0.5*r^2)
    */
    if (MAG_MAN(xrad) == 0)
    {
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
    }
    else
    {
        mag_t quadratic, comp_err;
        mp_limb_t A_sin, A_cos, A_exp;

        /* Bound computed error. */
        if (error != 0)
        {
            mag_init(comp_err); /* no need to free */
            mag_set_ui_2exp_si(comp_err, error, -wprounded);
        }

        /* Bound quadratic term for propagated error: 0.5*r^2 */
        mag_init(quadratic); /* no need to free */
        mag_nonzero_fast_mul(quadratic, xrad_copy, xrad_copy);
        MAG_EXP(quadratic) -= 1;

        /* Bound linear term for propagated error: cos(x)*r, sin(x)*r. */
        /* Note: we could have used the computed values, but then we would
           need to incorporate the computed error which would be slightly
           messier, and we would also need extra cases when only computing
           one of the functions. */
        /* Note: the bounds on cos(x) and sin(x) are assumed to be nonzero. */
        A_cos = arb_sin_cos_tab1[2 * p1_tab1][ARB_SIN_COS_TAB1_LIMBS - 1];
        A_sin = arb_sin_cos_tab1[2 * p1_tab1 + 1][ARB_SIN_COS_TAB1_LIMBS - 1];
        /* Adding 2 ulps (here ulp = 2^-8) gives an upper bound.
           The truncated table entry underestimates the sine or
           cosine of x by at most 1 ulp, and the top bits of x
           underestimate x by at most 1 ulp. */
        A_sin = (A_sin >> (FLINT_BITS - ARB_SIN_COS_TAB1_BITS)) + 2;
        A_cos = (A_cos >> (FLINT_BITS - ARB_SIN_COS_TAB1_BITS)) + 2;
        A_exp = -ARB_SIN_COS_TAB1_BITS;
        if (swapsincos)
        {
            mp_limb_t tt = A_sin;
            A_sin = A_cos;
            A_cos = tt;
        }
        /* Only use the top few bits */
        A_sin *= ((MAG_MAN(xrad_copy) >> (MAG_BITS - 16)) + LIMB_ONE);
        A_cos *= ((MAG_MAN(xrad_copy) >> (MAG_BITS - 16)) + LIMB_ONE);
        A_exp -= 16;
        A_exp += radexp;

        if (want_sin)
        {
            mag_set_ui_2exp_si(arb_radref(zsin), A_sin, A_exp);
            mag_nonzero_fast_add(arb_radref(zsin), arb_radref(zsin), quadratic);

            /* The propagated error is certainly at most r */
            if (mag_nonzero_fast_cmp(arb_radref(zsin), xrad_copy) > 0)
                mag_fast_set(arb_radref(zsin), xrad_copy);

            /* The propagated error is certainly at most 2 */
            if (MAG_EXP(arb_radref(zsin)) >= 2)
            {
                MAG_MAN(arb_radref(zsin)) = MAG_ONE_HALF;
                MAG_EXP(arb_radref(zsin)) = 2;
            }

            /* Add the computation error. */
            if (error != 0)
                mag_nonzero_fast_add(arb_radref(zsin), arb_radref(zsin), comp_err);
        }

        /* The same as above. */
        if (want_cos)
        {
            mag_set_ui_2exp_si(arb_radref(zcos), A_cos, A_exp);
            mag_nonzero_fast_add(arb_radref(zcos), arb_radref(zcos), quadratic);

            if (mag_nonzero_fast_cmp(arb_radref(zcos), xrad_copy) > 0)
                mag_fast_set(arb_radref(zcos), xrad_copy);

            if (MAG_EXP(arb_radref(zcos)) >= 2)
            {
                MAG_MAN(arb_radref(zcos)) = MAG_ONE_HALF;
                MAG_EXP(arb_radref(zcos)) = 2;
            }

            if (error != 0)
                mag_nonzero_fast_add(arb_radref(zcos), arb_radref(zcos), comp_err);
        }
    }

    /* Set the midpoints. */
    if (want_sin)
    {
        inexact = _arf_set_mpn_fixed(arb_midref(zsin), sinptr,
            wn, wn, sinnegative, prec, ARB_RND);
        if (inexact)
            arf_mag_add_ulp(arb_radref(zsin),
                arb_radref(zsin), arb_midref(zsin), prec);
    }

    if (want_cos)
    {
        inexact = _arf_set_mpn_fixed(arb_midref(zcos), cosptr,
            wn, wn, cosnegative, prec, ARB_RND);
        if (inexact)
            arf_mag_add_ulp(arb_radref(zcos),
                arb_radref(zcos), arb_midref(zcos), prec);
    }

    TMP_END;
}

void
arb_sin_cos(arb_t s, arb_t c, const arb_t x, slong prec)
{
    arb_sin_cos_fast(s, c, arb_midref(x), arb_radref(x), prec);
}

void
arb_sin(arb_t s, const arb_t x, slong prec)
{
    arb_sin_cos_fast(s, NULL, arb_midref(x), arb_radref(x), prec);
}

void
arb_cos(arb_t c, const arb_t x, slong prec)
{
    arb_sin_cos_fast(NULL, c, arb_midref(x), arb_radref(x), prec);
}

