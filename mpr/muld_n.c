#include "mpr.h"

/*
Computes z = x * y, rounding down (truncating) with correct rounding.
Returns adjustment (0 or 1) to be subtracted from the exponent.
*/

int
_mpr_muld_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    int shift;
    _MPR_MULD_N(shift, z, x, y, n);
    return shift;
}
