#include "mprb.h"

void
mprb_add(mprb_t z, const mprb_t x, const mprb_t y)
{
    mp_limb_t xrad, yrad, t, u;
    mp_srcptr xptr, yptr;
    long size, xexp, yexp, shift, carry, exp;

    if (x->sign != y->sign)
    {
        printf("mprb_add: different signs not supported\n");
        abort();
    }

    size = FLINT_MIN(x->size, y->size);
    size = FLINT_MIN(size, z->alloc);
    xrad = (size == x->size) ? x->rad : 2UL;
    yrad = (size == y->size) ? y->rad : 2UL;
    xptr = x->d + x->size - size;
    yptr = y->d + y->size - size;

    xexp = x->exp;
    yexp = y->exp;

    if (xexp >= yexp)
    {
        shift = xexp - yexp;
        carry = _mpr_addd_n(z->d, xptr, yptr, shift, size);
        exp = xexp;
    }
    else
    {
        shift = yexp - xexp;
        carry = _mpr_addd_n(z->d, yptr, xptr, shift, size);
        t = xrad;
        xrad = yrad;
        yrad = t;
        exp = yexp;
    }

    xrad = (xrad >> carry) + carry;
    yrad = (yrad >> carry) + carry;

    if (shift >= FLINT_BITS)
    {
        add_ssaaaa(u, t, 0, xrad, 0, 3);
    }
    else if (shift == 0)
    {
        add_ssaaaa(u, t, 0, xrad, 0, yrad);
        add_ssaaaa(u, t, u, t, 0, carry);
    }
    else
    {
        add_ssaaaa(u, t, 0, xrad, 0, (yrad >> shift) + 3);
    }

    if (u == 0)
    {
        z->rad = t;
        z->size = size;
        z->sign = x->sign;
        z->exp = exp + carry;
    }
    else
    {
        if (size == 1)
        {
            z->rad = 1UL << (FLINT_BITS - 1);
            z->d[0] = 0;
            z->size = 1;
            z->sign = x->sign;
            z->exp = exp + carry + 3;
        }
        else
        {
            z->rad = 3UL;
            mpn_copyi(z->d, z->d + 1, size - 1);
            z->size = size - 1;
            z->sign = x->sign;
            z->exp = exp + carry;
        }
    }
}
