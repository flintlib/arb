#include "mprb.h"

/* todo: handle zero */
void
mprb_set_mpfr(mprb_t x, const mpfr_t v)
{
    long xprec, vprec, prec;

    xprec = x->mid.alloc;
    vprec = _MPR_BITS_TO_LIMBS(v->_mpfr_prec);

    x->mid.exp = _mpr_set_mpfr(x->mid.d, v, xprec);
    x->mid.sign = (v->_mpfr_sign == 1) ? MPRB_SIGN_PLUS : MPRB_SIGN_MINUS;
    x->mid.size = x->mid.alloc;

    if (xprec >= vprec)
        ufloat_zero(&x->rad);
    else
        ufloat_set_2exp(&x->rad, mprb_ulp_exp(&x));
}
