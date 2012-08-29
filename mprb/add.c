#include "mprb.h"

void
mprb_add(mprb_t z, const mprb_t x, const mprb_t y)
{
    long error;

    ufloat_add(&z->rad, &x->rad, &y->rad);

    error = mpr_add_using_mpfr(&z->mid, &x->mid, &y->mid, z->bits, MPFR_RNDZ);

    if (error != LONG_MIN)
        mprb_add_rad_2exp(z, error);
}

/*
void
mprb_add(mprb_t z, const mprb_t x, const mprb_t y)
{
    mp_srcptr xptr, yptr;
    mp_ptr zptr;
    long zsize, xsize, ysize, xexp, yexp, shift;

    if (x->mid.sign != y->mid.sign)
        abort();

    xptr = x->mid.d;
    yptr = y->mid.d;
    zptr = z->mid.d;

    xsize = x->mid.size;
    ysize = y->mid.size;
    zsize = z->mid.alloc;

    xexp = x->mid.exp;
    yexp = y->mid.exp;

    ufloat_add(&z->rad, &x->rad, &y->rad);

    z->mid.sign = x->mid.sign;

    if (xsize == ysize && xsize == zsize)
    {
        if (xexp >= yexp)
        {
            shift = _mpr_addd_n(zptr, xptr, yptr, xexp - yexp, xsize);
            z->mid.exp = x->mid.exp + shift;
        }
        else
        {
            shift = _mpr_addd_n(zptr, yptr, xptr, yexp - xexp, xsize);
            z->mid.exp = y->mid.exp + shift;
        }
    }
    else
    {
        if (xexp >= yexp)
        {
            shift = _mpr_addd_using_mpfr(zptr, zsize, xptr, xsize, yptr, ysize, xexp - yexp);
            z->mid.exp = x->mid.exp + shift;
        }
        else
        {
            shift = _mpr_addd_using_mpfr(zptr, zsize, yptr, ysize, xptr, xsize, yexp - xexp);
            z->mid.exp = y->mid.exp + shift;
        }
    }

    z->mid.size = zsize;

    mprb_add_rad_2exp(z, mprb_ulp_exp(z));
}
*/