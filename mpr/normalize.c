#include "mpr.h"

/*
Normalization: given that the input mantissa is an arbitrary array of limbs,
shift and trim it so that the top bit of the top limb is set and
so that there are no trailing zero limbs, or if all limbs are zero,
set the input to the canonical representation of zero.

*/

/* remove any trailing zero limbs */
static __inline__ void
_mpr_trim(mpr_t x)
{
    long i = 0;

    while ((i < x->size) && (x->d[i] == 0UL)) i++;

    if (i != 0)
    {
        mpn_copyi(x->d, x->d + i, x->size - i);
        x->size = x->size - i;
    }
}

void
mpr_normalize(mpr_t x)
{
    long shift, b, size;
    mp_limb_t c;

    size = x->size;

    /* XXX: exponent if zero? */
    if (size == 0)
        return;

    for (shift = 0; shift < size; shift++)
    {
        c = x->d[size - 1 - shift];

        if (c != 0)
        {
            /* top bit is set */
            if (c & (1UL << (FLINT_BITS - 1)))
            {
                if (shift != 0)
                {
                    mpn_copyi(x->d, x->d + shift, size - shift);
                    x->size -= shift;
                }
            }
            else
            {
                count_leading_zeros(b, c);

                mpn_lshift(x->d, x->d, size, b);
                mpn_copyi(x->d, x->d + shift, size - shift);
                x->size -= shift;
                x->exp -= b;
            }

            break;
        }
    }

    _mpr_trim(x);

    /* XXX: exponent if zero? */
}
