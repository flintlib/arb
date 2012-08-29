#include "mpr.h"

int
mpr_is_normalized(mpr_t x)
{
    /* XXX: should look at the exponent */
    if (x->size == 0)
        return 1;

    /* there should be no trailing zero limbs */
    if (x->d[0] == 0UL)
        return 0;

    /* top bit should be set */
    if (!(x->d[x->size - 1] & (1UL << (FLINT_BITS - 1))))
        return 0;

    return 1;
}
