#include "mpr.h"

long
mpr_mul_using_mpfr(mpr_ptr z, mpr_srcptr x, mpr_srcptr y, mp_bitcnt_t prec, mpfr_rnd_t rnd)
{
    mpfr_t zs, xs, ys;
    long error;

    if (mpr_is_zero(x) || mpr_is_zero(y))
    {
        mpr_zero(z);
        return LONG_MIN;
    }

    prec = FLINT_MIN(prec, (x->size + y->size) * FLINT_BITS);

    mpr_fit_bits(z, prec);

    mpr_get_mpfr_wrapper(xs, x);
    mpr_get_mpfr_wrapper(ys, y);
    mpr_get_mpfr_wrapper(zs, z);

    zs->_mpfr_prec = prec;

    if (mpfr_mul(zs, xs, ys, rnd) != 0)
    {
        error = zs->_mpfr_exp - prec;
    }
    else
    {
        error = LONG_MIN;
    }

    z->size = _MPR_BITS_TO_LIMBS(prec);
    z->sign = zs->_mpfr_sign;
    z->exp = zs->_mpfr_exp;

    mpr_normalize(z);

    return error;
}

long
mpr_mul_ui_using_mpfr(mpr_ptr z, mpr_srcptr x, ulong y, mp_bitcnt_t prec, mpfr_rnd_t rnd)
{
    mpfr_t zs, xs;
    long error;

    if (mpr_is_zero(x) || y == 0)
    {
        mpr_zero(z);
        return LONG_MIN;
    }

    prec = FLINT_MIN(prec, (x->size + 1) * FLINT_BITS);

    mpr_fit_bits(z, prec);

    mpr_get_mpfr_wrapper(xs, x);
    mpr_get_mpfr_wrapper(zs, z);

    zs->_mpfr_prec = prec;

    if (mpfr_mul_ui(zs, xs, y, rnd) != 0)
    {
        error = zs->_mpfr_exp - prec;
    }
    else
    {
        error = LONG_MIN;
    }

    z->size = _MPR_BITS_TO_LIMBS(prec);
    z->sign = zs->_mpfr_sign;
    z->exp = zs->_mpfr_exp;

    mpr_normalize(z);

    return error;
}

long
mpr_mul_si_using_mpfr(mpr_ptr z, mpr_srcptr x, long y, mp_bitcnt_t prec, mpfr_rnd_t rnd)
{
    mpfr_t zs, xs;
    long error;

    if (mpr_is_zero(x) || y == 0)
    {
        mpr_zero(z);
        return LONG_MIN;
    }

    prec = FLINT_MIN(prec, (x->size + 1) * FLINT_BITS);

    mpr_fit_bits(z, prec);

    mpr_get_mpfr_wrapper(xs, x);
    mpr_get_mpfr_wrapper(zs, z);

    zs->_mpfr_prec = prec;

    if (mpfr_mul_si(zs, xs, y, rnd) != 0)
    {
        error = zs->_mpfr_exp - prec;
    }
    else
    {
        error = LONG_MIN;
    }

    z->size = _MPR_BITS_TO_LIMBS(prec);
    z->sign = zs->_mpfr_sign;
    z->exp = zs->_mpfr_exp;

    mpr_normalize(z);

    return error;
}
