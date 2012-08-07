#include "mprb.h"

void
mprb_mul(mprb_t z, const mprb_t x, const mprb_t y)
{
    ufloat_t a, b;
    mp_srcptr xptr, yptr;
    mp_ptr zptr;
    long zsize, xsize, ysize, exp, shift;

    xptr = x->mid.d;
    yptr = y->mid.d;
    zptr = z->mid.d;

    xsize = x->mid.size;
    ysize = y->mid.size;
    zsize = z->mid.size;

    _mpr_get_ufloat(a, xptr, xsize, x->mid.exp);
    _mpr_get_ufloat(b, yptr, ysize, y->mid.exp);

    /* error propagation: x*rad(y) + y*rad(x) + rad(x)*rad(y) */
    ufloat_mul(a, a, &y->rad);
    ufloat_addmul(a, b, &x->rad);
    ufloat_add_2exp(a, a, x->rad.exp + y->rad.exp);

    shift = _mpr_muld_using_mpfr(zptr, zsize, xptr, xsize, yptr, ysize);

    /* add exponents */
    exp = x->mid.exp + y->mid.exp - shift;

    /* rounding error */
    ufloat_add_2exp(a, a, exp - z->mid.size * FLINT_BITS);

    z->rad = *a;
    z->mid.exp = exp;
    z->mid.sign = x->mid.sign ^ y->mid.sign;
}
