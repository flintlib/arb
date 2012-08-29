#include "mpr.h"

/*
TODO:
* avoid temporary allocation
* make faster when ysize == 1

*/

#define TOP_BIT_SET(c) ((c) >> (FLINT_BITS-1))
#define MAX_STACK_ALLOC 32
#define MUL_BASECASE_CUTOFF 10

void
mpr_mul_exact(mpr_ptr z, mpr_srcptr x, mpr_srcptr y)
{
    long i, xsize, ysize, zsize;
    mp_limb_t tmp_stack[MAX_STACK_ALLOC];
    mp_ptr tmp;

    xsize = x->size;
    ysize = y->size;

    if (xsize < ysize)
    {
        mpr_srcptr t;
        long u;
        t = x; x = y; y = t;
        u = xsize;
        xsize = ysize;
        ysize = u;
    }

    if (y->size == 0)
    {
        mpr_zero(z);
        return;
    }

    zsize = xsize + ysize;

    if (zsize <= MAX_STACK_ALLOC)
    {
        tmp = tmp_stack;
    }
    else
    {
        tmp = malloc(sizeof(mp_limb_t) * zsize);
    }

    if (xsize < MUL_BASECASE_CUTOFF)
        mpn_mul_basecase(tmp, x->d, xsize, y->d, ysize);
    else
        mpn_mul(tmp, x->d, xsize, y->d, ysize);

    /* discard trailing zero limbs */
    i = 0;
    while (tmp[i] == 0)
        i++;

    /* remove limbs */
    if (TOP_BIT_SET(tmp[zsize - 1]))
    {
        mpr_fit_limbs(z, zsize - i);
        mpn_copyi(z->d, tmp + i, zsize - i);
        z->exp = x->exp + y->exp;
    }
    else
    {
        /* we might get one more zero limb after the shift */
        i += !(tmp[i] << 1);

        mpr_fit_limbs(z, zsize - i);
        mpn_lshift(z->d, tmp + i, zsize - i, 1);

        if (i > 0)
            z->d[0] |= (tmp[i - 1] >> (FLINT_BITS - 1));

        z->exp = x->exp + y->exp - 1;
    }

    if (zsize > MAX_STACK_ALLOC)
    {
        free(tmp);
    }

    z->size = zsize - i;
    z->sign = x->sign * y->sign;
}
