#include "mpr.h"


long
mpr_sub_using_mpfr(mpr_ptr z, mpr_srcptr x, mpr_srcptr y, mp_bitcnt_t prec, mpfr_rnd_t rnd)
{
    mpfr_t zs, xs, ys;
    long error;

    prec = FLINT_MIN(prec,
        1 + FLINT_MAX(x->exp, y->exp) - FLINT_MIN(x->exp - x->size * FLINT_BITS,
                                y->exp - y->size * FLINT_BITS));

    mpr_fit_bits(z, prec);

    mpr_get_mpfr_wrapper(xs, x);
    mpr_get_mpfr_wrapper(ys, y);
    mpr_get_mpfr_wrapper(zs, z);

    zs->_mpfr_prec = prec;

    if (mpfr_sub(zs, xs, ys, rnd) != 0)
    {
        error = zs->_mpfr_exp - prec;
    }
    else
    {
        error = LONG_MIN;
    }

    if (mpfr_zero_p(zs))
    {
        mpr_zero(z);
        return error;
    }

    z->size = _MPR_BITS_TO_LIMBS(prec);
    z->sign = zs->_mpfr_sign;
    z->exp = zs->_mpfr_exp;

    mpr_normalize(z);

    return error;
}
