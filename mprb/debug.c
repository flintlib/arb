#include "mprb.h"

void
mprb_debug(const mprb_t x)
{
    long e;
    mpz_t t;

    mpz_init(t);
    e = mprb_get_mid_mpz_2exp(t, x);

    printf("N: [exp=%ld] ", x->exp);
    mpn_debug(x->d, x->size);
    gmp_printf("Z: {mid=%Zx, exp=%ld, rad=%lu, size=%ld, alloc=%ld}\n",
        t, e, x->rad, (long) x->size, (long) x->alloc);

    mpz_clear(t);
}
