#include "mprb.h"

/* XXX: signs! */
void
mprb_get_mid_mpfr(mpfr_t x, const mprb_t v, mpfr_rnd_t rnd)
{
    if ((v->size == 1) && (v->d[0] == 0))
        mpfr_set_ui(x, 0, MPFR_RNDD);
    else
        _mpr_get_mpfr_signed(x, v->d, v->exp, v->size,
            (v->sign == MPRB_SIGN_PLUS) ? 1 : -1, rnd);
}
