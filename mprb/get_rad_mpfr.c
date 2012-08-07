#include "mprb.h"

void
mprb_get_rad_mpfr(mpfr_t r, const mprb_t x)
{
    ufloat_get_mpfr(r, &x->rad);
}
