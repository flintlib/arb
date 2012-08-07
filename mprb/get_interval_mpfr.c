#include "mprb.h"

void
mprb_get_interval_mpfr(mpfr_t a, mpfr_t b, const mprb_t x)
{
    mpfr_t r;
    mpfr_init2(r, FLINT_BITS);

    mprb_get_mid_mpfr(a, x, MPFR_RNDD);
    mprb_get_mid_mpfr(b, x, MPFR_RNDU);

    mprb_get_rad_mpfr(r, x);

    mpfr_sub(a, a, r, MPFR_RNDD);
    mpfr_add(b, b, r, MPFR_RNDU);

    mpfr_clear(r);
}
