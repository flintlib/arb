#include "mpr.h"

#define _MPR_MUL_SHIFT(shift, d, dn, s, sn) \
    do { \
        int i; \
        shift = !((s)[(sn)-1] >> (FLINT_BITS-1)); \
        for (i = (dn) - 1; i >= 0; i--) \
            (d)[i] = ((s)[i+(sn)-(dn)] << shift) \
                   | (((s)[i+(sn)-(dn)-1] >> (FLINT_BITS-1)) & shift); \
    } while (0)

#define _MPR_MULD_N_CONST(shift, z, x, y, n) \
    do { \
        mp_limb_t __tmp[2*n]; \
        MPN_MUL_N(__tmp, x, y, n); \
        _MPR_MUL_SHIFT(shift, z, n, __tmp, 2*n); \
    } while (0)

#define _MPR_MULD_N(shift, z, x, y, n) \
    do { \
        if (n == 1) \
            _MPR_MULD_N_CONST(shift, z, x, y, 1); \
        else if (n == 2) \
            _MPR_MULD_N_CONST(shift, z, x, y, 2); \
        else if (n <= 15) \
        { \
            mp_limb_t __tmp[30]; \
            mpn_mul_basecase(__tmp, x, n, y, n); \
            _MPR_MUL_SHIFT(shift, z, n, __tmp, 2*n); \
        } \
        else \
        { \
            mp_ptr __tmp = malloc(2 * n * sizeof(mp_limb_t)); \
            mpn_mul_n(__tmp, x, y, n); \
            _MPR_MUL_SHIFT(shift, z, n, __tmp, 2*n); \
            free(__tmp); \
        } \
    } while (0)


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
