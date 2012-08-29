#include "mpr.h"

void
_mpr_randtest(mp_ptr x, flint_rand_t state, mp_size_t n)
{
    __mpz_struct z;
    mp_limb_t t;

    _flint_rand_init_gmp(state);

    z._mp_size = n;
    z._mp_alloc = n;
    z._mp_d = x;

    t = n_randlimb(state);

    if (t & 1)
        mpz_urandomb(&z, state->gmp_state, n * FLINT_BITS);
    else
        mpz_rrandomb(&z, state->gmp_state, n * FLINT_BITS);

    x[n - 1] |= (1UL << (FLINT_BITS - 1));
}

void
mpr_randtest(mpr_t x, flint_rand_t state, long bits)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_randtest(t, state, bits);
    mpr_set_fmpz_2exp(x, t, n_randint(state, 20) - 10);
    fmpz_clear(t);
    mpr_normalize(x);
}
