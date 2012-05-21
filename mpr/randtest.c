#include "mpr.h"

void
_mpr_randtest(mp_ptr x, flint_rand_t state, mp_size_t n)
{
    __mpz_struct z;

    _flint_rand_init_gmp(state);

    z._mp_size = n;
    z._mp_alloc = n;
    z._mp_d = x;

if (0)    mpz_urandomb(&z, state->gmp_state, n * FLINT_BITS);

    mpn_random2(x, n);

    x[n - 1] |= (1UL << (FLINT_BITS - 1));
}
