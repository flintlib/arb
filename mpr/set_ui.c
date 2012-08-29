#include "mpr.h"

void mpr_set_ui(mpr_t y, ulong x)
{
    if (x == 0)
    {
        mpr_zero(y);
    }
    else
    {
        long c;

        mpr_fit_limbs(y, 1);

        count_leading_zeros(c, x);

        y->d[0] = x << c;
        y->exp = FLINT_BITS - c;
        y->size = 1;
        y->sign = MPR_SIGN_PLUS;
    }
}
