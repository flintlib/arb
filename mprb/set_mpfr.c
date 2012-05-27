#include "mprb.h"

/* todo: handle zero */
void
mprb_set_mpfr(mprb_t x, const mpfr_t v)
{
    long xprec, vprec, prec;

    xprec = x->alloc;
    vprec = _MPR_BITS_TO_LIMBS(v->_mpfr_prec);

    x->exp = _mpr_set_mpfr(x->d, v, xprec);
    x->rad = (xprec >= vprec) ? 0UL : 1UL;
    x->sign = (v->_mpfr_sign == 1) ? MPRB_SIGN_PLUS : MPRB_SIGN_MINUS;
    x->size = x->alloc;
}
