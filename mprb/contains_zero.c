#include "mprb.h"

int
mprb_contains_zero(const mprb_t x)
{
    mp_limb_t midm, radm;
    long i;

    /* centered on zero */
    if (x->mid.d[x->mid.size - 1] == 0)
        return 1;

    /* midpoint is larger */
    if (x->mid.exp > x->rad.exp)
        return 0;

    /* radius is larger */
    if (x->mid.exp < x->rad.exp)
        return 1;

    /* same exponent; compare mantissas */
    midm = x->mid.d[x->mid.size - 1];
    radm = x->rad.man << (FLINT_BITS - UFLOAT_PREC);

    /* midpoint is larger */
    if (midm > radm)
        return 0;

    /* radius is larger */
    if (midm < radm)
        return 1;

    /* midpoint is larger */
    for (i = x->mid.size - 2; i >= 0; i--)
    {
        if (x->mid.d[i] != 0)
            return 0;
    }

    /* equal */
    return 1;
}
