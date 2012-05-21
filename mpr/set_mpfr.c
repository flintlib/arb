#include "mpr.h"

long
_mpr_set_mpfr(mp_ptr y, const mpfr_t x, mp_size_t n, mpfr_rnd_t rnd)
{
    long xlimbs, exp, m;

    if (!mpfr_regular_p(x) || mpfr_sgn(x) != 1)
    {
        printf("_mpr_set_mpfr: input must be a positive regular number\n");
        abort();
    }

    xlimbs = _MPR_BITS_TO_LIMBS(mpfr_get_prec(x));

    if (n >= xlimbs || rnd == MPFR_RNDD || rnd == MPFR_RNDZ)
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
        mpfr_set(t, x, rnd);
        mpn_copyi(y, t->_mpfr_d, n);
        exp = mpfr_get_exp(t);
        mpfr_clear(t);
    }

    return exp;
}
