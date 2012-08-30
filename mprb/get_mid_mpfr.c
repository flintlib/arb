#include "mprb.h"

/* XXX: signs! */
void
mprb_get_mid_mpfr(mpfr_t x, const mprb_t v, mpfr_rnd_t rnd)
{
    mpr_get_mpfr(x, v, rnd);
}
