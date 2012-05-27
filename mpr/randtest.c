#include "mpr.h"

void
_mpr_randtest(mp_ptr x, flint_rand_t state, mp_size_t n)
{
    __mpz_struct z;
    mp_limb_t t;
    long i;

    _flint_rand_init_gmp(state);

    z._mp_size = n;
    z._mp_alloc = n;
    z._mp_d = x;

    t = n_randlimb(state);

    if (t & 1)
        mpz_urandomb(&z, state->gmp_state, n * FLINT_BITS);
    else
        mpz_rrandomb(&z, state->gmp_state, n * FLINT_BITS);

    for (i = z._mp_size; i < n; i++)
        x[i] = 0;

    x[n - 1] |= (1UL << (FLINT_BITS - 1));
}
