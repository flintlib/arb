#include "mprb.h"

/* todo: handle 0 */
long
mprb_get_mid_mpz_2exp(mpz_t a, const mprb_t x)
{
    if ((x->size == 1) && (x->d[0] == 0))
    {
        mpz_set_ui(a, 0);
        return 0;
    }
    else
    {
        mpz_realloc2(a, x->size * FLINT_BITS);
        mpn_copyi(a->_mp_d, x->d, x->size);
        a->_mp_size = (x->sign == MPRB_SIGN_PLUS) ? x->size : -(x->size);
        return x->exp - x->size * FLINT_BITS;
    }
}
