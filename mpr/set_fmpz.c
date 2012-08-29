#include "mpr.h"

void mpr_set_fmpz(mpr_t y, const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
    {
        mpr_set_si(y, *x);
    }
    else
    {
        mpr_set_mpz(y, COEFF_TO_PTR(*x));
    }
}
