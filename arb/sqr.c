/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

#define ARB_TRIM_BITS (MAG_BITS / 2)

#define SQR_MPFR_MIN_LIMBS 40
#define SQR_MPFR_MAX_LIMBS 10000

static __inline__ int
_arf_mul_1x1(arf_t res, slong * exp_shift, mp_limb_t xtop, mp_limb_t ytop, int negative, slong prec)
{
    mp_limb_t hi, lo;
    int inexact;

    /* The result will have 1 or 2 limbs */
    ARF_DEMOTE(res);

    umul_ppmm(hi, lo, xtop, ytop);

    /* Align top bit. Should this be branchless? */
    /* todo: use comparison with LIMB_TOP instead of shift? */
    if (!(hi >> (FLINT_BITS - 1)))
    {
        hi = (hi << 1) | (lo >> (FLINT_BITS - 1));
        lo = (lo << 1);
        *exp_shift = -1;
    }
    else
    {
        *exp_shift = 0;
    }

    if (prec <= FLINT_BITS)
    {
        mp_limb_t hh;
        hh = hi & (LIMB_ONES << (FLINT_BITS - prec));
        inexact = (hh != hi) || (lo != 0);
        hi = hh;
        ARF_NOPTR_D(res)[0] = hi;
        ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, negative);
    }
    else
    {
        if (prec <= 2 * FLINT_BITS)
        {
            mp_limb_t ll;
            ll = lo & (LIMB_ONES << (2 * FLINT_BITS - prec));
            inexact = (ll != lo);
            lo = ll;
        }
        else
        {
            inexact = 0;
        }

        if (lo != 0)
        {
            ARF_NOPTR_D(res)[1] = hi;
            ARF_NOPTR_D(res)[0] = lo;
            ARF_XSIZE(res) = ARF_MAKE_XSIZE(2, negative);
        }
        else
        {
            ARF_NOPTR_D(res)[0] = hi;
            ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, negative);
        }
    }

    return inexact;
}

static __inline__ int
_arf_mul_2x1(arf_t res, slong * exp_shift, mp_limb_t xtop, mp_limb_t xlo, mp_limb_t ytop, int negative, slong prec)
{
    mp_limb_t tmp[3];
    mp_limb_t hi, lo;
    int inexact;

    nn_mul_2x1(tmp[2], tmp[1], tmp[0], xtop, xlo, ytop);

    if (prec <= 2 * FLINT_BITS)
    {
        ARF_DEMOTE(res);

        /* First align and truncate to 2 limbs */
        if ((tmp[2] >> (FLINT_BITS - 1)) == 0)
        {
            tmp[2] = (tmp[2] << 1) | (tmp[1] >> (FLINT_BITS - 1));
            tmp[1] = (tmp[1] << 1) | (tmp[0] >> (FLINT_BITS - 1));
            inexact = ((tmp[0] << 1) != 0);
            *exp_shift = -1;
        }
        else
        {
            inexact = (tmp[0] != 0);
            *exp_shift = 0;
        }

        hi = tmp[2];
        lo = tmp[1];

        if (prec <= FLINT_BITS)
        {
            mp_limb_t hh;
            hh = hi & (LIMB_ONES << (FLINT_BITS - prec));
            inexact = inexact || (lo != 0) || (hh != hi);
            hi = hh;
            ARF_NOPTR_D(res)[0] = hi;
            ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, negative);
        }
        else
        {
            mp_limb_t ll;
            ll = lo & (LIMB_ONES << (2 * FLINT_BITS - prec));
            inexact = inexact || (ll != lo);
            lo = ll;

            if (lo != 0)
            {
                ARF_NOPTR_D(res)[1] = hi;
                ARF_NOPTR_D(res)[0] = lo;
                ARF_XSIZE(res) = ARF_MAKE_XSIZE(2, negative);
            }
            else
            {
                ARF_NOPTR_D(res)[0] = hi;
                ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, negative);
            }
        }
    }
    else
    {
        inexact = _arf_set_round_mpn(res, exp_shift, tmp, 3, negative, prec, ARF_RND_DOWN);
    }

    return inexact;
}

static __inline__ int
_arf_mul_2x2(arf_t res, slong * exp_shift, mp_limb_t xtop, mp_limb_t xlo, mp_limb_t ytop, mp_limb_t ylo, int negative, slong prec)
{
    mp_limb_t tmp[4];
    mp_limb_t hi, lo;
    int inexact;

    nn_mul_2x2(tmp[3], tmp[2], tmp[1], tmp[0], xtop, xlo, ytop, ylo);

    if (prec <= 2 * FLINT_BITS)
    {
        ARF_DEMOTE(res);

        /* First align and truncate to 2 limbs */
        if ((tmp[3] >> (FLINT_BITS - 1)) == 0)
        {
            tmp[3] = (tmp[3] << 1) | (tmp[2] >> (FLINT_BITS - 1));
            tmp[2] = (tmp[2] << 1) | (tmp[1] >> (FLINT_BITS - 1));
            inexact = (tmp[0] != 0) || ((tmp[1] << 1) != 0);
            *exp_shift = -1;
        }
        else
        {
            inexact = (tmp[0] != 0) || (tmp[1] != 0);
            *exp_shift = 0;
        }

        hi = tmp[3];
        lo = tmp[2];

        if (prec <= FLINT_BITS)
        {
            mp_limb_t hh;
            hh = hi & (LIMB_ONES << (FLINT_BITS - prec));
            inexact = inexact || (lo != 0) || (hh != hi);
            hi = hh;
            ARF_NOPTR_D(res)[0] = hi;
            ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, negative);
        }
        else
        {
            mp_limb_t ll;
            ll = lo & (LIMB_ONES << (2 * FLINT_BITS - prec));
            inexact = inexact || (ll != lo);
            lo = ll;

            if (lo != 0)
            {
                ARF_NOPTR_D(res)[1] = hi;
                ARF_NOPTR_D(res)[0] = lo;
                ARF_XSIZE(res) = ARF_MAKE_XSIZE(2, negative);
            }
            else
            {
                ARF_NOPTR_D(res)[0] = hi;
                ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, negative);
            }
        }
    }
    else
    {
        inexact = _arf_set_round_mpn(res, exp_shift, tmp, 4, negative, prec, ARF_RND_DOWN);
    }

    return inexact;
}

int
_arf_mul_via_mpfr(arf_ptr z, slong * exp_shift,
    mp_srcptr xptr, mp_size_t xn, mp_srcptr yptr, mp_size_t yn,
    int negative, slong prec, arf_rnd_t rnd)
{
    mp_size_t zn, val;
    mp_ptr tmp, zptr;
    mpfr_t xf, yf, zf;
    int ret;
    ARF_MUL_TMP_DECL

    prec = FLINT_MIN((xn + yn) * FLINT_BITS, prec);
    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    ARF_MUL_TMP_ALLOC(tmp, zn)

    zf->_mpfr_d = tmp;
    zf->_mpfr_prec = prec;
    zf->_mpfr_sign = 1;
    zf->_mpfr_exp = 0;

    xf->_mpfr_d = (mp_ptr) xptr;
    xf->_mpfr_prec = xn * FLINT_BITS;
    xf->_mpfr_sign = negative ? -1 : 1;
    xf->_mpfr_exp = 0;

    if (xptr == yptr && xn == yn)
    {
        ret = mpfr_sqr(zf, xf, arf_rnd_to_mpfr(rnd));
    }
    else
    {
        yf->_mpfr_d = (mp_ptr) yptr;
        yf->_mpfr_prec = yn * FLINT_BITS;
        yf->_mpfr_sign = 0;
        yf->_mpfr_exp = 0;

        ret = mpfr_mul(zf, xf, yf, arf_rnd_to_mpfr(rnd));
    }

    ret = (ret != 0);
    *exp_shift = zf->_mpfr_exp;

    val = 0;
    while (tmp[val] == 0)
        val++;

    ARF_GET_MPN_WRITE(zptr, zn - val, z);
    flint_mpn_copyi(zptr, tmp + val, zn - val);
    ARF_XSIZE(z) |= negative;

    ARF_MUL_TMP_FREE(tmp, zn)

    return ret;
}

/* Is it worth having a special version of this for small exponents? */
static void
_arb_sqr_no_prec(arb_t res, const arb_t x)
{
    arb_get_mag(arb_radref(res), x);
    mag_mul(arb_radref(res), arb_radref(res), arb_radref(res));
    arf_zero(arb_midref(res));
}

/* todo: 3x3 sqr when prec <= 192? */

int
_arf_sqr_mpn(arf_t res, slong * exp_shift, mp_srcptr x, mp_size_t xn, slong prec)
{
    int inexact;

    if (xn <= 2)
    {
        if (xn == 1)
            inexact = _arf_mul_1x1(res, exp_shift, x[0], x[0], 0, prec);
        else
            inexact = _arf_mul_2x2(res, exp_shift, x[1], x[0], x[1], x[0], 0, prec);
    }
    /* MPFR implements mulhigh which speeds up squaring in a certain range. */
    else if (xn >= SQR_MPFR_MIN_LIMBS &&
             xn <= SQR_MPFR_MAX_LIMBS &&
             prec < xn * (7 * FLINT_BITS / 4))   /* Computing much less than the full product. */
    {
        inexact = _arf_mul_via_mpfr(res, exp_shift, x, xn, x, xn, 0, prec, ARF_RND_DOWN);
    }
    else
    {
        mp_ptr tmp;
        ARF_MUL_TMP_DECL
        ARF_MUL_TMP_ALLOC(tmp, 2 * xn)
        mpn_sqr(tmp, x, xn);
        inexact = _arf_set_round_mpn(res, exp_shift, tmp, 2 * xn, 0, prec, ARF_RND_DOWN);
        ARF_MUL_TMP_FREE(tmp, 2 * xn)
    }

    return inexact;
}

/* Slow path: infinities and/or large exponents. */
void
arb_sqr_fallback(arb_t res, const arb_t x, slong prec)
{
    if (!arb_is_finite(x))
    {
        /*
        (nan +/- c)^2    ->  nan +/- inf
        (inf +/- c)^2    ->  inf +/- 0
        (-inf +/- c)^2   ->  inf +/- 0
        (c +/- inf)^2    ->  0 +/- inf
        */

        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else if (mag_is_finite(arb_radref(x)))
            arb_pos_inf(res);
        else
            arb_zero_pm_inf(res);
    }
    else if (arf_is_zero(arb_midref(x)))
    {
        arf_zero(arb_midref(res));
        mag_mul(arb_radref(res), arb_radref(x), arb_radref(x));
    }
    else
    {
        mp_size_t xn, xuse;
        mp_srcptr xptr;
        int inexact;
        slong accuracy, accuracy_limbs, exp_shift;

        ARF_GET_MPN_READONLY(xptr, xn, arb_midref(x));

        if (mag_is_zero(arb_radref(x)))
        {
            accuracy = MAG_MAX_LAGOM_EXP;
        }
        else
        {
            accuracy = _fmpz_sub_small(ARF_EXPREF(arb_midref(x)), MAG_EXPREF(arb_radref(x)));
            accuracy = FLINT_MIN(accuracy, MAG_MAX_LAGOM_EXP);
        }
        accuracy = FLINT_MIN(accuracy, prec) + ARB_TRIM_BITS;

        prec = FLINT_MIN(prec, accuracy);
        accuracy_limbs = (accuracy + FLINT_BITS - 1) / FLINT_BITS;
        xuse = FLINT_MIN(xn, accuracy_limbs);

        if (prec >= 2)
        {
            mag_t xrad;
            mag_init(xrad);

            /* Add to xrad the error due to truncating the midpoint. */
            if (xuse != xn)
                arf_mag_add_ulp(xrad, arb_radref(x), arb_midref(x), FLINT_BITS * xuse);
            else
                mag_set(xrad, arb_radref(x));

            arf_get_mag(arb_radref(res), arb_midref(x));

            /* (x+xrad)^2 = x^2 + 2x xrad + xrad^2 */
            mag_mul(arb_radref(res), arb_radref(res), xrad);
            mag_mul_2exp_si(arb_radref(res), arb_radref(res), 1);
            mag_addmul(arb_radref(res), xrad, xrad);

            inexact = _arf_sqr_mpn(arb_midref(res), &exp_shift, xptr + xn - xuse, xuse, prec);

            _fmpz_add2_fast(ARF_EXPREF(arb_midref(res)),
                            ARF_EXPREF(arb_midref(x)),
                            ARF_EXPREF(arb_midref(x)), exp_shift);

            if (inexact)
                arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec);

            mag_clear(xrad);
        }
        else
        {
            _arb_sqr_no_prec(res, x);
        }
    }
}

void
arb_sqr(arb_t res, const arb_t x, slong prec)
{
    mp_size_t xn, xuse;
    mp_srcptr xptr;
    mp_limb_t xtop;
    unsigned int xrad;
    slong xexp, zexp, xradexp, accuracy, accuracy_limbs;
    slong exp_shift;
    int inexact;

    if (!ARB_IS_LAGOM(x))
    {
        arb_sqr_fallback(res, x, prec);
        return;
    }

    /* Demote any bignum exponents in the output. This does nothing if the
       exponents are small, so aliasing is safe. */
    _fmpz_demote(ARF_EXPREF(arb_midref(res)));
    _fmpz_demote(MAG_EXPREF(arb_radref(res)));

    xrad = MAG_MAN(arb_radref(x));
    xn = ARF_SIZE(arb_midref(x));

    /* The midpoint is zero. */
    if (xn == 0)
    {
        /* Zero the output midpoint. */
        ARF_DEMOTE(arb_midref(res));
        ARF_XSIZE(arb_midref(res)) = 0;
        ARF_EXP(arb_midref(res)) = ARF_EXP_ZERO;

        /* Square the radius. */
        if (xrad == 0)
        {
            MAG_MAN(arb_radref(res)) = 0;
            MAG_EXP(arb_radref(res)) = 0;
        }
        else
        {
            xradexp = MAG_EXP(arb_radref(x));
            MAG_MUL(xrad, xradexp, xrad, xradexp, xrad, xradexp);
            MAG_MAN(arb_radref(res)) = xrad;
            MAG_EXP(arb_radref(res)) = xradexp;
        }
        return;
    }

    xexp = ARF_EXP(arb_midref(x));
    zexp = 2 * xexp;

    /* Fast path for low-precision input. */
    if (xn <= 2)
    {
        if (xrad == 0)  /* Exact input (no error propagation). */
        {
            /* Multiply the midpoints. */
            if (xn == 1)
            {
                xtop = ARF_NOPTR_D(arb_midref(x))[0];
                inexact = _arf_mul_1x1(arb_midref(res), &exp_shift, xtop, xtop, 0, prec);
            }
            else
            {
                mp_limb_t xlo;
                xtop = ARF_NOPTR_D(arb_midref(x))[1];
                xlo = ARF_NOPTR_D(arb_midref(x))[0];
                inexact = _arf_mul_2x2(arb_midref(res), &exp_shift, xtop, xlo, xtop, xlo, 0, prec);
            }

            zexp += exp_shift;
            ARF_EXP(arb_midref(res)) = zexp;

            /* Write the radius. */
            if (inexact)
            {
                MAG_SET_2EXP(MAG_MAN(arb_radref(res)), MAG_EXP(arb_radref(res)), zexp - prec);
            }
            else
            {
                MAG_MAN(arb_radref(res)) = 0;
                MAG_EXP(arb_radref(res)) = 0;
            }
        }
        else   /* We need to take the input error into account. */
        {
            /* Calculate lower bound for the output accuracy and reduce the
               output precision accordingly. */
            xradexp = MAG_EXP(arb_radref(x));
            accuracy = xexp - xradexp + ARB_TRIM_BITS;
            prec = FLINT_MIN(prec, accuracy);

            if (prec >= 2)
            {
                if (xn == 1)
                {
                    xtop = ARF_NOPTR_D(arb_midref(x))[0];
                    inexact = _arf_mul_1x1(arb_midref(res), &exp_shift, xtop, xtop, 0, prec);
                }
                else
                {
                    /* Note: even if we have reduced the needed precision to 1 limb,
                       and xn == 2, it is cheap enough to just do the 2x2 multiply
                       instead of fiddling around with the error bounds. */
                    mp_limb_t xlo;
                    xtop = ARF_NOPTR_D(arb_midref(x))[1];
                    xlo = ARF_NOPTR_D(arb_midref(x))[0];
                    inexact = _arf_mul_2x2(arb_midref(res), &exp_shift, xtop, xlo, xtop, xlo, 0, prec);
                }

                zexp += exp_shift;
                ARF_EXP(arb_midref(res)) = zexp;

                /* Propagated error: (x +/- xrad)^2 -> x^2 +/- (2 |x| xrad + xrad^2). */
                MAG_SET_ARF_TOP_LIMB(xtop, xexp, xtop, xexp);
                MAG_MULADDMUL(xtop, xexp, xtop, xexp + 1, xrad, xradexp, xrad, xradexp, xrad, xradexp);

                /* Add rounding error. */
                if (inexact)
                    MAG_ADD_2EXP(xtop, xexp, xtop, xexp, zexp - prec);

                /* Finally, write the radius. */
                MAG_MAN(arb_radref(res)) = xtop;
                MAG_EXP(arb_radref(res)) = xexp;
            }
            else
            {
                _arb_sqr_no_prec(res, x);
            }
        }
    }
    else
    {
        xptr = ARF_PTR_D(arb_midref(x));
        xtop = xptr[xn - 1];
        xradexp = MAG_EXP(arb_radref(x));

        /* Calculate lower bound for the output accuracy and reduce the
           output precision accordingly. */
        if (xrad != 0)
            accuracy = xexp - xradexp;
        else
            accuracy = MAG_MAX_LAGOM_EXP;
        accuracy = FLINT_MIN(accuracy, prec) + ARB_TRIM_BITS;
        prec = FLINT_MIN(prec, accuracy);

        /* Throw away insignificant input limbs */
        accuracy_limbs = (accuracy + FLINT_BITS - 1) / FLINT_BITS;
        xuse = FLINT_MIN(xn, accuracy_limbs);

        if (prec >= 2)
        {
            inexact = _arf_sqr_mpn(arb_midref(res), &exp_shift, xptr + xn - xuse, xuse, prec);
            zexp += exp_shift;
            ARF_EXP(arb_midref(res)) = zexp;

            /* Add to xrad the error due to truncating the midpoint. */
            if (xuse != xn)
            {
                if (xrad == 0)
                    MAG_SET_2EXP(xrad, xradexp, xexp - FLINT_BITS * xuse);
                else
                    MAG_ADD_2EXP(xrad, xradexp, xrad, xradexp, xexp - FLINT_BITS * xuse);
            }

            if (xrad == 0)
            {
                if (inexact)
                {
                    MAG_SET_2EXP(xtop, xexp, zexp - prec);
                    MAG_MAN(arb_radref(res)) = xtop;
                    MAG_EXP(arb_radref(res)) = xexp;
                }
                else
                {
                    MAG_MAN(arb_radref(res)) = 0;
                    MAG_EXP(arb_radref(res)) = 0;
                }
            }
            else
            {
                MAG_SET_ARF_TOP_LIMB(xtop, xexp, xtop, xexp);
                MAG_MULADDMUL(xtop, xexp, xtop, xexp + 1, xrad, xradexp, xrad, xradexp, xrad, xradexp);

                /* Add rounding error. */
                if (inexact)
                {
                    MAG_ADD_2EXP(xtop, xexp, xtop, xexp, zexp - prec);
                }

                /* Finally, write the radius. */
                MAG_MAN(arb_radref(res)) = xtop;
                MAG_EXP(arb_radref(res)) = xexp;
            }
        }
        else
        {
            _arb_sqr_no_prec(res, x);
        }
    }
}

