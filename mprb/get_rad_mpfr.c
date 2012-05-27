#include "mprb.h"

void
mprb_get_rad_mpfr(mpfr_t r, const mprb_t x)
{
    mpfr_set_ui_2exp(r, x->rad, x->exp - FLINT_BITS * x->size, MPFR_RNDU);
}
