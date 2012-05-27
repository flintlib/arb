#include "mprb.h"

void
mprb_get_interval_mpfr(mpfr_t a, mpfr_t b, const mprb_t x)
{
    mpfr_t r;
    mpfr_init2(r, FLINT_BITS);

    mprb_get_mid_mpfr(a, x, MPFR_RNDD);
    mprb_get_mid_mpfr(b, x, MPFR_RNDU);

    mprb_get_rad_mpfr(r, x);

#if 0
    printf("GETTING INTRVLA\n");
    mpfr_printf("a = %Rg [xsign = %ld]\n", a, x->sign);
    mpfr_printf("b = %Rg\n", b);
    mpfr_printf("r = %Rg\n", r);
#endif

    mpfr_sub(a, a, r, MPFR_RNDD);
    mpfr_add(b, b, r, MPFR_RNDU);

#if 0
    mpfr_printf("a2 = %Rg\n", a);
    mpfr_printf("b2 = %Rg\n", b);
#endif

    mpfr_clear(r);
}
