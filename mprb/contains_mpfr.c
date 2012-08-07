#include "mprb.h"

int
mprb_contains_mpfr(const mprb_t x, const mpfr_t f)
{
    mpfr_t a, b;
    int result;

    /* XXX: this is no longer correct if the radius is larger */

    mpfr_init2(a, FLINT_BITS * x->mid.size + 1000);
    mpfr_init2(b, FLINT_BITS * x->mid.size + 1000);

    mprb_get_interval_mpfr(a, b, x);

    result = (mpfr_cmp(a, f) <= 0) && (mpfr_cmp(f, b) <= 0);

    mpfr_clear(a);
    mpfr_clear(b);

    return result;
}
