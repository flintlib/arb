#include "mpr.h"

void mpr_set_mpz(mpr_t y, const mpz_t x)
{
    long c, size;

    size = x->_mp_size;

    if (size == 0)
    {
        mpr_zero(y);
    }
    else
    {
        size = FLINT_ABS(size);

        mpr_fit_limbs(y, size);

        count_leading_zeros(c, x->_mp_d[size - 1]);

        if (c == 0)
            mpn_copyi(y->d, x->_mp_d, size);
        else
            mpn_lshift(y->d, x->_mp_d, size, c);

        y->size = size;
        y->exp = size * FLINT_BITS - c;
        y->sign = (x->_mp_size < 0) ? MPR_SIGN_MINUS : MPR_SIGN_PLUS;

        mpr_normalize(y);
    }
}
