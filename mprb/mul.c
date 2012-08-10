#include "mprb.h"

void
mprb_mul(mprb_t z, const mprb_t x, const mprb_t y)
{
    ufloat_t a, b;
    mp_srcptr xptr, yptr;
    mp_ptr zptr;
    long zsize, xsize, ysize, shift;

    xptr = x->mid.d;
    yptr = y->mid.d;
    zptr = z->mid.d;

    xsize = x->mid.size;
    ysize = y->mid.size;
    zsize = z->mid.alloc;


    /* error propagation: x*rad(y) + y*rad(x) + rad(x)*rad(y) */
#if 0
    ufloat_set_2exp(a, a, x->mid.exp + y->rad.exp);
    ufloat_add_2exp(a, a, y->mid.exp + x->rad.exp);
    ufloat_add_2exp(a, a, x->rad.exp + y->rad.exp);
#else
    mpr_get_ufloat(a, &x->mid);
    mpr_get_ufloat(b, &y->mid);
    ufloat_mul(a, a, &y->rad);
    ufloat_addmul(a, b, &x->rad);
    ufloat_add_2exp(a, a, x->rad.exp + y->rad.exp);
#endif

    if (xsize == ysize && xsize == zsize && xsize <= 8)
        shift = _mpr_muld_n(zptr, xptr, yptr, xsize);
    else
        shift = _mpr_muld_using_mpfr(zptr, zsize, xptr, xsize, yptr, ysize);

    /* add exponents */
    z->mid.exp = x->mid.exp + y->mid.exp - shift;

    /* rounding error */
    ufloat_add_2exp(a, a, mprb_ulp_exp(z));

    z->rad = *a;
    z->mid.size = zsize;
    z->mid.sign = x->mid.sign ^ y->mid.sign;
}


#if 0

void mprb_mul(mprb_t z, const mprb_t x, const mprb_t y)
{
    ufloat_t a, error;
    long xsize, ysize, zsize, zmax, rounding_exp;

    if (mprb_is_special(x) || mprb_is_special(y))
    {
        abort();
    }

    xsize = mprb_mid(x)->size;
    ysize = mprb_mid(y)->size;
    zsize = xsize + ysize;
    zmax = mprb_mid(z)->max;

    if (mprb_is_exact(x))
    {
        if (mprb_is_exact(y))
        {
            mpr_mul_exact(mprb_mid(z), mprb_mid(x), mprb_mid(y));
            ufloat_zero(mprb_rad(z));
            return;
        }
        else
        {
            /* propagated error: x*rad(y) */
            zsize = FLINT_MIN(ysize, zmax);
            mpr_get_ufloat(error, mprb_mid(x));
            ufloat_mul(error, error, mprb_rad(y));
        }
    }
    else
    {
        if (mprb_is_exact(y))
        {
            /* propagated error: y*rad(x) */
            zsize = FLINT_MIN(xsize, zmax);
            mpr_get_ufloat(error, mprb_mid(y));
            ufloat_mul(error, error, mprb_rad(x));
        }
        else
        {
            zsize = FLINT_MIN(xsize, ysize);
            zsize = FLINT_MIN(zsize, zmax);

            /* propagated error: x*rad(y) + y*rad(x) + rad(x)*rad(y) */
            mpr_get_ufloat(a, mprb_mid(x));
            mpr_get_ufloat(b, mprb_mid(y));

            ufloat_mul(a, a, mprb_rad(y));
            ufloat_addmul(a, b, mprb_rad(x));

            /* rad(x)*rad(y) is usually tiny, so approximate it */
            ufloat_add_2exp(error, error, mprb_rad_exp(x) + mprb_rad_exp(y));
        }
    }

    rounding_exp = mpr_mul_trunc(mprb_mid(z), mprb_mid(x), mprb_mid(y), zsize);

    if (rounding_exp == LONG_MIN)
        ufloat_set(mprb_rad(z), error);
    else
        ufloat_add_2exp(mprb_rad(z), error, err_exp);

    /* todo: adjust precision to accuracy */

}

#endif
