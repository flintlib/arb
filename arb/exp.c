/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

#define TMP_ALLOC_LIMBS(__n) TMP_ALLOC((__n) * sizeof(mp_limb_t))
#define MAGLIM(prec) FLINT_MAX(128, 2 * (prec))

void
arb_exp_arf_huge(arb_t z, const arf_t x, slong mag, slong prec, int minus_one)
{
    arb_t ln2, t, u;
    fmpz_t q;
    slong wp;

    arb_init(ln2);
    arb_init(t);
    arb_init(u);
    fmpz_init(q);

    wp = prec + mag + 10;

    arb_const_log2(ln2, wp);
    arb_set_arf(t, x);
    arb_div(u, t, ln2, wp);
    arf_get_fmpz(q, arb_midref(u), ARF_RND_DOWN);
    arb_submul_fmpz(t, ln2, q, wp);

    arb_exp(z, t, prec);
    arb_mul_2exp_fmpz(z, z, q);

    if (minus_one)
        arb_sub_ui(z, z, 1, prec);

    arb_clear(ln2);
    arb_clear(t);
    arb_clear(u);
    fmpz_clear(q);
}

/* |x| >= 2^expbound */
static void
arb_exp_arf_overflow(arb_t z, slong expbound, int negative, int minus_one, slong prec)
{
    if (!negative)
    {
        arf_zero(arb_midref(z));
        mag_inf(arb_radref(z));
    }
    else
    {
        /* x <= -2^expbound   ==>   0 < exp(x) <= 2^(-2^expbound) */
        fmpz_t t;
        fmpz_init(t);

        fmpz_set_si(t, -1);
        fmpz_mul_2exp(t, t, expbound);

        arf_one(arb_midref(z));
        mag_one(arb_radref(z));
        arb_mul_2exp_fmpz(z, z, t);

        if (minus_one)
            arb_sub_ui(z, z, 1, prec);

        fmpz_clear(t);
    }
}

static void
arb_exp_arf_fallback(arb_t z, const arf_t x, slong mag, slong prec, int minus_one)
{
    /* reduce by log(2) if needed, but avoid computing log(2) unnecessarily at
       extremely high precision */
    if (mag > 64 || (mag > 8 && prec < 1000000))
        arb_exp_arf_huge(z, x, mag, prec, minus_one);
    else if (prec < 19000)
        arb_exp_arf_rs_generic(z, x, prec, minus_one);
    else
        arb_exp_arf_bb(z, x, prec, minus_one);
}

void
arb_exp_arf(arb_t z, const arf_t x, slong prec, int minus_one, slong maglim)
{
    if (arf_is_special(x))
    {
        if (minus_one)
        {
            if (arf_is_zero(x))
                arb_zero(z);
            else if (arf_is_pos_inf(x))
                arb_pos_inf(z);
            else if (arf_is_neg_inf(x))
                arb_set_si(z, -1);
            else
                arb_indeterminate(z);
        }
        else
        {
            if (arf_is_zero(x))
                arb_one(z);
            else if (arf_is_pos_inf(x))
                arb_pos_inf(z);
            else if (arf_is_neg_inf(x))
                arb_zero(z);
            else
                arb_indeterminate(z);
        }
    }
    else if (COEFF_IS_MPZ(ARF_EXP(x)))
    {
        if (fmpz_sgn(ARF_EXPREF(x)) > 0)
        {
            /* huge input */
            arb_exp_arf_overflow(z, maglim, ARF_SGNBIT(x), minus_one, prec);
        }
        else
        {
            /* |exp(x) - (1 + x)| <= |x^2| */
            fmpz_t t;
            int inexact;
            fmpz_init(t);
            fmpz_mul_2exp(t, ARF_EXPREF(x), 1);
            inexact = arf_add_ui(arb_midref(z), x, minus_one ? 0 : 1, prec, ARB_RND);
            mag_one(arb_radref(z));
            mag_mul_2exp_fmpz(arb_radref(z), arb_radref(z), t);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
            fmpz_clear(t);
        }
    }
    else
    {
        slong exp, wp, wn, N, r, wprounded, finaln;
        fmpz_t n;
        mp_ptr tmp, w, t, u, finalvalue;
        mp_limb_t p1, q1bits, p2, q2bits, error, error2;
        int negative, inexact;
        TMP_INIT;

        exp = ARF_EXP(x);
        negative = ARF_SGNBIT(x);

        /* handle tiny input */
        /* |exp(x) - 1| <= 2|x| */
        if (!minus_one && exp < -prec - 4)
        {
            arf_one(arb_midref(z));
            mag_set_ui_2exp_si(arb_radref(z), 1, exp + 1);
            return;
        }
        /* |exp(x) - (1 + x)| <= |x^2| */
        else if (exp < (minus_one ? -prec - 4 : -(prec / 2) - 4))
        {
            inexact = arf_add_ui(arb_midref(z), x, minus_one ? 0 : 1, prec, ARB_RND);
            mag_set_ui_2exp_si(arb_radref(z), 1, 2 * exp);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
            return;
        }

        /* handle huge input */
        if (exp > maglim)
        {
            arb_exp_arf_overflow(z, maglim, negative, minus_one, prec);
            return;
        }

        /* Absolute working precision (NOT rounded to a limb multiple) */
        wp = prec + 8;
        if (minus_one && exp <= 0)
            wp += (-exp);
        /* Number of limbs */
        wn = (wp + FLINT_BITS - 1) / FLINT_BITS;
        /* Precision rounded to a number of bits */
        wprounded = FLINT_BITS * wn;
        /* Don't be close to the boundary (to allow adding adding the
           Taylor series truncation error without overflow) */
        wp = FLINT_MAX(wp, wprounded - (FLINT_BITS - 4));

        /* Too high precision to use table -- use generic algorithm */
        if (wp > ARB_EXP_TAB2_PREC)
        {
            arb_exp_arf_fallback(z, x, exp, prec, minus_one);
            return;
        }

        TMP_START;

        tmp = TMP_ALLOC_LIMBS(4 * wn + 3);
        w = tmp;        /* requires wn+1 limbs */
        t = w + wn + 1; /* requires wn+1 limbs */
        u = t + wn + 1; /* requires 2wn+1 limbs */

        /* reduce modulo log(2) */
        fmpz_init(n);
        if (_arb_get_mpn_fixed_mod_log2(w, n, &error, x, wn) == 0)
        {
            /* may run out of precision for log(2) */
            arb_exp_arf_fallback(z, x, exp, prec, minus_one);
            fmpz_clear(n);
            TMP_END;
            return;
        }

        /* err(w) translates to a propagated error bounded by
           err(w) * exp'(x) < err(w) * exp(1) < err(w) * 3 */
        error *= 3;

        /* Table-based argument reduction (1 or 2 steps) */
        if (wp <= ARB_EXP_TAB1_PREC)
        {
            q1bits = ARB_EXP_TAB1_BITS;
            q2bits = 0;

            p1 = w[wn-1] >> (FLINT_BITS - q1bits);
            w[wn-1] -= (p1 << (FLINT_BITS - q1bits));
            p2 = 0;
        }
        else
        {
            q1bits = ARB_EXP_TAB21_BITS;
            q2bits = ARB_EXP_TAB21_BITS + ARB_EXP_TAB22_BITS;

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

        if (N < 60)
        {
            /* Evaluate Taylor series */
            _arb_exp_taylor_rs(t, &error2, w, wn, N);
            /* Taylor series evaluation error */
            error += error2;
            /* Taylor series truncation error */
            error += UWORD(1) << (wprounded-wp);
        }
        else  /* Compute cosh(a) from sinh(a) using a square root. */
        {
            /* the summation for sinh is actually done to (2N-1)! */
            N = (N + 1) / 2;

            /* Evaluate Taylor series for sinh */
            _arb_sin_cos_taylor_rs(t, u, &error2, w, wn, N, 1, 0);
            error += error2;
            error += UWORD(1) << (wprounded-wp);

            /* 1 + sinh^2, with wn + 1 limbs */
            mpn_sqr(u, t, wn);
            u[2 * wn] = 1;

            /* cosh, with wn + 1 limbs */
            mpn_sqrtrem(w, u, u, 2 * wn + 1);

            /* exp = sinh + cosh */
            t[wn] = w[wn] + mpn_add_n(t, t, w, wn);

            /* Error for cosh */
            /* When converting sinh to cosh, the error for cosh must be
               smaller than the error for sinh; but we also get 1 ulp
               extra error from the square root. */
            error2 = error + 1;

            /* Error for sinh + cosh */
            error += error2;
        }

        if (wp <= ARB_EXP_TAB1_PREC)
        {
            if (p1 == 0)
            {
                finalvalue = t;
                finaln = wn + 1;
            }
            else
            {
                /* Divide by 2 to get |t| <= 1 (todo: check this?) */
                mpn_rshift(t, t, wn + 1, 1);
                error = (error >> 1) + 2;

                mpn_mul_n(u, t, arb_exp_tab1[p1] + ARB_EXP_TAB1_LIMBS - wn, wn);

                /* (t + err1 * ulp) * (u + err2 * ulp) + 1ulp = t*u +
                   (err1*u + err2*t + t*u*ulp + 1) * ulp
                   note |u| <= 1, |t| <= 1 */
                error += 4;
                finalvalue = u + wn;

                finaln = wn;

                /* we have effectively divided by 2^2 -- todo use inline function */
                fmpz_add_ui(n, n, 2);
            }
        }
        else
        {
            if (p1 == 0 && p2 == 0)
            {
                finalvalue = t;
                finaln = wn + 1;
            }
            else
            {
                /* Divide by 2 to get |t| <= 1 (todo: check this?) */
                mpn_rshift(t, t, wn + 1, 1);
                error = (error >> 1) + 2;

                mpn_mul_n(u, arb_exp_tab21[p1] + ARB_EXP_TAB2_LIMBS - wn,
                             arb_exp_tab22[p2] + ARB_EXP_TAB2_LIMBS - wn, wn);

                /* error of w <= 4 ulp */
                flint_mpn_copyi(w, u + wn, wn);  /* todo: avoid with better alloc */

                mpn_mul_n(u, t, w, wn);

                /* (t + err1 * ulp) * (w + 4 * ulp) + 1ulp = t*u +
                   (err1*w + 4*t + t*w*ulp + 1) * ulp
                   note |w| <= 1, |t| <= 1 */
                error += 6;

                finalvalue = u + wn;
                finaln = wn;

                /* we have effectively divided by 2^3 -- todo use inline function */
                fmpz_add_ui(n, n, 3);
            }
        }

        /* The accumulated arithmetic error */
        mag_set_ui_2exp_si(arb_radref(z), error, -wprounded);

        /* Set the midpoint */
        if (!minus_one)
        {
            inexact = _arf_set_mpn_fixed(arb_midref(z), finalvalue, finaln, wn, 0, prec, ARB_RND);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
        }
        else
        {
            _arf_set_mpn_fixed(arb_midref(z), finalvalue, finaln, wn, 0, finaln * FLINT_BITS, ARB_RND);
        }

        arb_mul_2exp_fmpz(z, z, n);

        if (minus_one)
            arb_sub_ui(z, z, 1, prec);

        TMP_END;
        fmpz_clear(n);
    }
}

/* todo: min prec by MAG_BITS everywhere? */
void
arb_exp_wide(arb_t res, const arb_t x, slong prec, slong maglim)
{
    mag_t t, u;

    mag_init(t);
    mag_init(u);

    if (arf_cmpabs_2exp_si(arb_midref(x), 20) < 0
        && mag_cmp_2exp_si(arb_radref(x), 20) < 0)
    {
        if (arf_is_zero(arb_midref(x)))
        {
            if (mag_cmp_2exp_si(arb_radref(x), -10) < 0)
            {
                mag_expm1(arb_radref(res), arb_radref(x));
                arf_one(arb_midref(res));
            }
            else
            {
                mag_expinv_lower(t, arb_radref(x));
                mag_exp(u, arb_radref(x));
                arb_set_interval_mag(res, t, u, prec);
            }
        }
        else if (arb_contains_zero(x))
        {
            arf_get_mag_lower(t, arb_midref(x));
            mag_sub(t, arb_radref(x), t);

            arf_get_mag(u, arb_midref(x));
            mag_add(u, arb_radref(x), u);

            if (arf_sgn(arb_midref(x)) > 0)
            {
                mag_expinv_lower(t, t);
                mag_exp(u, u);
                arb_set_interval_mag(res, t, u, prec);
            }
            else
            {
                mag_expinv_lower(u, u);
                mag_exp(t, t);
                arb_set_interval_mag(res, u, t, prec);
            }
        }
        else if (arf_sgn(arb_midref(x)) < 0)
        {
            arb_get_mag(t, x);
            arb_get_mag_lower(u, x);
            mag_expinv_lower(t, t);
            mag_expinv(u, u);
            arb_set_interval_mag(res, t, u, prec);
        }
        else
        {
            arb_get_mag_lower(t, x);
            arb_get_mag(u, x);
            mag_exp_lower(t, t);
            mag_exp(u, u);
            arb_set_interval_mag(res, t, u, prec);
        }
    }
    else
    {
        /* use arb_exp_arf for accurate argument reduction */
        arf_t q;
        arf_init(q);
        arf_set_mag(q, arb_radref(x));
        arf_add(q, q, arb_midref(x), MAG_BITS, ARF_RND_CEIL);
        arb_exp_arf(res, q, FLINT_MIN(prec, MAG_BITS), 0, maglim);
        arb_get_mag(arb_radref(res), res);
        arf_zero(arb_midref(res));
        arf_clear(q);
    }

    mag_clear(t);
    mag_clear(u);
}

void arb_exp(arb_t res, const arb_t x, slong prec)
{
    slong maglim = MAGLIM(prec);

    if (mag_is_zero(arb_radref(x)))
    {
        arb_exp_arf(res, arb_midref(x), prec, 0, maglim);
    }
    else if (mag_is_inf(arb_radref(x)))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else
            arb_zero_pm_inf(res);
    }
    else if (arf_is_special(arb_midref(x)))
    {
        if (arf_is_zero(arb_midref(x)))
            arb_exp_wide(res, x, prec, maglim);
        else if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else  /* infinity +/- finite */
            arb_exp_arf(res, arb_midref(x), prec, 0, 1);
    }
    else  /* both finite, non-special */
    {
        slong acc, mexp, rexp;

        mexp = ARF_EXP(arb_midref(x));
        rexp = MAG_EXP(arb_radref(x));

        if (COEFF_IS_MPZ(rexp))
            rexp = (fmpz_sgn(ARF_EXPREF(arb_radref(x))) < 0) ? COEFF_MIN : COEFF_MAX;
        if (COEFF_IS_MPZ(mexp))
            mexp = (fmpz_sgn(MAG_EXPREF(arb_midref(x))) < 0) ? COEFF_MIN : COEFF_MAX;

        if (mexp < -prec && rexp < -prec)
        {
            arb_get_mag(arb_radref(res), x);
            mag_expm1(arb_radref(res), arb_radref(res));
            arf_one(arb_midref(res));
            return;
        }

        acc = -rexp;
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc < 20 && (rexp >= 0 || mexp <= 10))
        {
            /* may evaluate at endpoints */
            arb_exp_wide(res, x, prec, maglim);
        }
        else
        {
            /* exp(a+b) - exp(a) = exp(a) * (exp(b)-1) */
            mag_t t, u;

            mag_init_set(t, arb_radref(x));
            mag_init(u);

            arb_exp_arf(res, arb_midref(x), prec, 0, maglim);
            mag_expm1(t, t);
            arb_get_mag(u, res);
            mag_addmul(arb_radref(res), t, u);

            mag_clear(t);
            mag_clear(u);
        }
    }
}

void
arb_expm1(arb_t res, const arb_t x, slong prec)
{
    slong maglim = MAGLIM(prec);

    if (mag_is_zero(arb_radref(x)))
    {
        arb_exp_arf(res, arb_midref(x), prec, 1, maglim);
    }
    else if (mag_is_inf(arb_radref(x)))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else
            arb_zero_pm_inf(res);
    }
    else if (arf_is_special(arb_midref(x)))
    {
        if (arf_is_zero(arb_midref(x)))
        {
            if (mag_cmp_2exp_si(arb_radref(x), -10) < 0)
            {
                mag_expm1(arb_radref(res), arb_radref(x));
                arf_zero(arb_midref(res));
            }
            else
            {
                arb_exp_wide(res, x, prec, maglim);
                arb_sub_ui(res, res, 1, prec);
            }
        }
        else if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else  /* infinity +/- finite */
            arb_exp_arf(res, arb_midref(x), prec, 1, 1);
    }
    else  /* both finite, non-special */
    {
        if (arf_cmpabs_2exp_si(arb_midref(x), 3) < 0 &&
            mag_cmp_2exp_si(arb_radref(x), -3) < 0)
        {
            mag_t t, u, one;
            slong acc, mexp, rexp;

            mexp = ARF_EXP(arb_midref(x));
            rexp = MAG_EXP(arb_radref(x));

            if (COEFF_IS_MPZ(rexp))
                rexp = (fmpz_sgn(ARF_EXPREF(arb_radref(x))) < 0) ? COEFF_MIN : COEFF_MAX;
            if (COEFF_IS_MPZ(mexp))
                mexp = (fmpz_sgn(MAG_EXPREF(arb_midref(x))) < 0) ? COEFF_MIN : COEFF_MAX;

            acc = FLINT_MIN(mexp, 0) - rexp;
            acc = FLINT_MAX(acc, 0);
            acc = FLINT_MIN(acc, prec);
            prec = FLINT_MIN(prec, acc + MAG_BITS);
            prec = FLINT_MAX(prec, 2);

            /* [exp(a+b) - 1] - [exp(a) - 1] = exp(a) * (exp(b)-1) */
            mag_init_set(t, arb_radref(x));
            mag_init(u);
            mag_init(one);
            mag_one(one);

            if (arf_sgn(arb_midref(x)) >= 0)
            {
                arb_exp_arf(res, arb_midref(x), prec, 1, maglim);
                arb_get_mag(u, res);
                mag_add(u, u, one);
            }
            else
            {
                arb_exp_arf(res, arb_midref(x), prec, 1, maglim);
                arb_get_mag_lower(u, res);
                mag_sub(u, one, u);
            }

            mag_expm1(t, t);
            mag_addmul(arb_radref(res), t, u);

            mag_clear(t);
            mag_clear(u);
            mag_clear(one);
        }
        else
        {
            arb_exp(res, x, prec);
            arb_sub_ui(res, res, 1, prec);
        }
    }
}

void arb_exp_invexp(arb_t res, arb_t res2, const arb_t x, slong prec)
{
    slong maglim = MAGLIM(prec);

    if (arf_is_special(arb_midref(x)) || mag_is_special(arb_radref(x)))
    {
        /* [c +/- 0] */
        if (arf_is_finite(arb_midref(x)) && mag_is_zero(arb_radref(x)))
        {
            arb_exp_arf(res, arb_midref(x), prec, 0, maglim);
            arb_inv(res2, res, prec);
        }  /* [nan +/- ?] */
        else if (arf_is_nan(arb_midref(x)))
        {
            arb_indeterminate(res);
            arb_indeterminate(res2);
        }  /* [c +/- inf] */
        else if (mag_is_inf(arb_radref(x)))
        {
            arb_zero_pm_inf(res);
            arb_zero_pm_inf(res2);
        }  /* [+inf +/- c] */
        else if (arf_is_pos_inf(arb_midref(x)))
        {
            arb_pos_inf(res);
            arb_zero(res2);
        }  /* [-inf +/- c] */
        else if (arf_is_neg_inf(arb_midref(x)))
        {
            arb_zero(res);
            arb_pos_inf(res2);
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, x);
            arb_exp_wide(res, x, prec, maglim);
            arb_exp_wide(res2, t, prec, maglim);
            arb_clear(t);
        }
    }
    else  /* both finite, non-special */
    {
        slong acc, mexp, rexp;

        mexp = ARF_EXP(arb_midref(x));
        rexp = MAG_EXP(arb_radref(x));

        if (COEFF_IS_MPZ(rexp))
            rexp = (fmpz_sgn(ARF_EXPREF(arb_radref(x))) < 0) ? COEFF_MIN : COEFF_MAX;
        if (COEFF_IS_MPZ(mexp))
            mexp = (fmpz_sgn(MAG_EXPREF(arb_midref(x))) < 0) ? COEFF_MIN : COEFF_MAX;

        if (mexp < -prec && rexp < -prec)
        {
            arb_get_mag(arb_radref(res), x);
            mag_expm1(arb_radref(res), arb_radref(res));
            arf_one(arb_midref(res));
            arb_set(res2, res);
            return;
        }

        acc = -rexp;
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc < 20 && (rexp >= 0 || mexp <= 10))
        {
            /* may evaluate at endpoints */
            arb_t t;
            arb_init(t);
            arb_neg(t, x);
            arb_exp_wide(res, x, prec, maglim);
            arb_exp_wide(res2, t, prec, maglim);
            arb_clear(t);
        }
        else
        {
            /* exp(a+b) - exp(a) = exp(a) * (exp(b)-1) */
            mag_t t, u;

            mag_init_set(t, arb_radref(x));
            mag_init(u);

            arb_exp_arf(res, arb_midref(x), prec, 0, maglim);
            arb_inv(res2, res, prec);

            mag_expm1(t, t);

            arb_get_mag(u, res);
            mag_addmul(arb_radref(res), t, u);
            arb_get_mag(u, res2);
            mag_addmul(arb_radref(res2), t, u);

            mag_clear(t);
            mag_clear(u);
        }
    }
}
