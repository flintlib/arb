#include "mprb.h"

void
mprb_debug(const mprb_t x)
{
    long e;
    mpz_t t;

    mpz_init(t);
    e = mprb_get_mid_mpz_2exp(t, x);

/*
    printf("N: [exp=%ld] ", x->mid.exp);
    mpn_debug(x->mid.d, x->mid.size);
*/

    gmp_printf("{%Zd * 2^%ld +- %lu * 2^%ld, size=%ld, alloc=%ld}\n",
        t, e, x->rad.man, x->rad.exp - UFLOAT_PREC, (long) x->mid.size, (long) x->mid.alloc);

    mpz_clear(t);
}
