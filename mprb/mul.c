#include "mprb.h"

/* throw away to avoid overflow */
#define RAD_MUL_SHIFT 8

void
mprb_mul(mprb_t z, const mprb_t x, const mprb_t y)
{
    mp_limb_t xrad, yrad;
    mp_srcptr xptr, yptr;
    long size;

    size = FLINT_MIN(x->size, y->size);
    size = FLINT_MIN(size, z->alloc);
    xrad = (size == x->size) ? x->rad : 2UL;
    yrad = (size == y->size) ? y->rad : 2UL;
    xptr = x->d + x->size - size;
    yptr = y->d + y->size - size;

    if (size == 1)
    {
        mp_limb_t u1, u0, a2, a1, a0, b1, b0;

        umul_ppmm(a1, a0, xptr[0], yrad);
        umul_ppmm(b1, b0, yptr[0], xrad);
        add_sssaaaaaa(a2, a1, a0, 0, a1, a0, 0, b1, b0);
        umul_ppmm(b1, b0, xrad, yrad);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, b1, b0);

        umul_ppmm(u1, u0, xptr[0], yptr[0]);

        if ((u1 >> (FLINT_BITS - 1)) == 0)
        {
            u1 = (u1 << 1) | (u0 >> (FLINT_BITS - 1));
            u0 = (u0 << 1);
            a2 = (a2 << 1) | (a1 >> (FLINT_BITS - 1));
            a1 = (a1 << 1) | (a0 >> (FLINT_BITS - 1));
            a0 = (a0 << 1);
            z->exp = x->exp + y->exp - 1;
        }
        else
        {
            z->exp = x->exp + y->exp;
        }

        /* round up */
        add_ssaaaa(a2, a1, a2, a1, 0, (u0 != 0) + (a0 != 0));

        if (a2 == 0)
        {
            z->rad = a1;
            z->d[0] = u1;
        }
        else
        {
            z->rad = -1UL;
            z->d[0] = (1UL << (FLINT_BITS - 1));
            z->exp += 2;
        }

        z->size = 1;
        z->sign = x->sign ^ y->sign;
    }
    else
    {
        mp_limb_t t0, t1, u0, u1, rad;
        mp_limb_t xtop, ytop;
        int shift, need_limb_shift;

        /* upper bound for x*yrad + y*xrad */
        xtop = (xptr[size-1] >> RAD_MUL_SHIFT) + 1UL;
        ytop = (yptr[size-1] >> RAD_MUL_SHIFT) + 1UL;
        umul_ppmm(t1, t0, xtop, y->rad);
        umul_ppmm(u1, u0, ytop, x->rad);
        add_ssaaaa(t1, t0, t1, t0, u1, u0);

        /* scale back radius to the true magnitude */
        t0 = (t0 >> (FLINT_BITS - RAD_MUL_SHIFT)) | (t1 << RAD_MUL_SHIFT);
        t1 = (t1 >> (FLINT_BITS - RAD_MUL_SHIFT));

        /* xrad * yrad < 1 ulp
            + 1 ulp error in the main multiplication
            + 1 ulp error from the shift above */
        add_ssaaaa(t1, t0, t1, t0, 0, 3UL);

        shift = _mpr_muld_n(z->d, xptr, yptr, size);

        /* adjust radius */
        t1 = (t1 << shift) | ((t0 >> (FLINT_BITS - 1)) & shift);
        t0 = t0 << shift;

        /* remove one limb */
        if (t1 != 0)
        {
            need_limb_shift = 1;

            /* +1 ulp error for the rounding of t1 */
            /* +1 ulp error in the main product, which doubles
                if we perform a shift */
            rad = t1 + 2 + shift;
            mpn_copyi(z->d, z->d + 1, size - 1);
        }
        else
        {
            need_limb_shift = 0;
            rad = t0;
        }

        z->rad = rad;
        z->exp = x->exp + y->exp - shift;
        z->size = size - need_limb_shift;
        z->sign = x->sign ^ y->sign;
    }
}
