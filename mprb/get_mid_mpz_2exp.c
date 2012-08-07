#include "mprb.h"

/* todo: handle 0 */
long
mprb_get_mid_mpz_2exp(mpz_t a, const mprb_t x)
{
    if ((x->mid.size == 1) && (x->mid.d[0] == 0))
    {
        mpz_set_ui(a, 0);
        return 0;
    }
    else
    {
        mpz_realloc2(a, x->mid.size * FLINT_BITS);
        mpn_copyi(a->_mp_d, x->mid.d, x->mid.size);
        a->_mp_size = (x->mid.sign == MPRB_SIGN_PLUS) ? x->mid.size : -(x->mid.size);
        return x->mid.exp - x->mid.size * FLINT_BITS;
    }
}
