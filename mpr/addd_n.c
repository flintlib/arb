#include "mpr.h"

#define ADD_MAX_STACK_ALLOC 256

/*
Computes z = x + y / 2^shift, rounding down (truncating) with correct rounding.
Returns adjustment (0 or 1) to be added to the exponent.
*/

int
_mpr_addd_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_bitcnt_t shift, mp_size_t n)
{
    mp_limb_t cy;

    if (shift >= n * FLINT_BITS)
    {
        if (z != x)
            mpn_copyi(z, x, n);
        return 0;
    }

    if (n == 1)
    {
        mp_limb_t t;
        t = x[0] + (y[0] >> shift);
        cy = t < x[0];
        z[0] = (t >> cy) | (cy << (FLINT_BITS - 1));
    }
    else if (n == 2)
    {
        mp_limb_t s1, s0, y1, y0;

        y1 = y[1];
        y0 = y[0];

        if (shift < FLINT_BITS)
        {
            if (shift != 0)
            {
                y0 = (y0 >> shift) | (y1 << (FLINT_BITS - shift));
                y1 = y1 >> shift;
            }
        }
        else
        {
            y0 = y1 >> shift;
            y1 = 0;
        }

        add_sssaaaaaa(cy, s1, s0, 0, x[1], x[0], 0, y1, y0);

        z[1] = (s1 >> cy) | (cy << (FLINT_BITS - 1));
        z[0] = (s0 >> cy) | ((s1 & cy) << (FLINT_BITS - 1));
    }
    else
    {
        mp_size_t b, s, m;

        b = shift % FLINT_BITS;
        s = shift / FLINT_BITS;
        m = n - s;

        if (b == 0)
        {
            cy = mpn_add_n(z, x, y + s, m);
            if (s & cy)
                cy = mpn_add_1(z + m, z + m, 1, s);
        }
        else
        {
            if (m <= ADD_MAX_STACK_ALLOC)
            {
                mp_limb_t tmp[ADD_MAX_STACK_ALLOC];
                mpn_rshift(tmp, y + s, m, b);
                cy = mpn_add_n(z, x, tmp, m);
            }
            else
            {
                mp_ptr tmp = malloc(sizeof(mp_limb_t) * m);
                mpn_rshift(tmp, y + s, m, b);
                cy = mpn_add_n(z, x, tmp, m);
                free(tmp);
            }

            if (s & cy)
                cy = mpn_add_1(z + n - s, z + n - s, 1, s);
        }

        if (cy)
        {
            mpn_rshift(z, z, n, 1);
            z[n - 1] |= (1UL << (FLINT_BITS - 1));
        }
    }

    return cy;
}

