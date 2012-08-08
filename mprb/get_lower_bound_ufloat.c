#include "mprb.h"

/* sets u to a lower bound for mid(x)-rad(x), assuming that this is positive

XXX: implementation is wrong if rad is zero
*/

void
mprb_get_lower_bound_ufloat(ufloat_t u, const mprb_t x)
{
    __mpfr_struct zs, xs, ys;

    mp_limb_t r, s;

    xs._mpfr_d = (mp_ptr) x->mid.d;
    xs._mpfr_exp = x->mid.exp;
    xs._mpfr_prec = x->mid.size * FLINT_BITS;
    xs._mpfr_sign = 1;

    r = x->rad.man << (FLINT_BITS - UFLOAT_PREC);

    ys._mpfr_d = &r;
    ys._mpfr_exp = x->rad.exp;
    ys._mpfr_prec = UFLOAT_PREC;
    ys._mpfr_sign = 1;

    zs._mpfr_d = &s;
    zs._mpfr_exp = 0;
    zs._mpfr_prec = UFLOAT_PREC;
    zs._mpfr_sign = 1;

    mpfr_sub(&zs, &xs, &ys, MPFR_RNDD);

    u->man = s >> (FLINT_BITS - UFLOAT_PREC);
    u->exp = zs._mpfr_exp;
}
