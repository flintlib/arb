#include "mpr.h"

long
_mpr_set_mpfr(mp_ptr y, const mpfr_t x, mp_size_t n)
{
    long xlimbs, exp, m;

    if (!mpfr_regular_p(x))
    {
        printf("_mpr_set_mpfr: input must be a regular number\n");
        abort();
    }

    xlimbs = _MPR_BITS_TO_LIMBS(mpfr_get_prec(x));

    if (n >= xlimbs)
    {
        m = FLINT_MIN(xlimbs, n);
        mpn_copyi(y + n - m, x->_mpfr_d + m - xlimbs, m);
        mpn_zero(y, n - m);
        exp = mpfr_get_exp(x);
    }
    else
    {
        mpfr_t t;
        mpfr_init2(t, n * FLINT_BITS);
        mpfr_set(t, x, MPFR_RNDZ);
        mpn_copyi(y, t->_mpfr_d, n);
        exp = mpfr_get_exp(t);
        mpfr_clear(t);
    }

    return exp;
}

void mpr_set_mpfr(mpr_t y, const mpfr_t x)
{
    if (mpfr_zero_p(x))
    {
        mpr_zero(y);
    }
    else
    {
        long n = _MPR_BITS_TO_LIMBS(mpfr_get_prec(x));

        mpr_fit_limbs(y, n);

        y->exp = _mpr_set_mpfr(y->d, x, n);
        y->size = n;
        y->sign = x->_mpfr_sign;

        mpr_normalize(y);
    }
}
