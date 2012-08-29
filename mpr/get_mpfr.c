#include "mpr.h"

void
_mpr_get_mpfr(mpfr_t y, mp_srcptr x, long exp, mp_size_t n, mpfr_rnd_t rnd)
{
    __mpz_struct z;

    z._mp_size = n;
    z._mp_alloc = n;
    z._mp_d = (mp_ptr) x;

    mpfr_set_z_2exp(y, &z, exp - n * FLINT_BITS, rnd);
}

void
_mpr_get_mpfr_signed(mpfr_t y, mp_srcptr x, long exp,
    mp_size_t n, int sign, mpfr_rnd_t rnd)
{
    __mpz_struct z;

    z._mp_size = (sign > 0) ? n : -n;
    z._mp_alloc = n;
    z._mp_d = (mp_ptr) x;

    mpfr_set_z_2exp(y, &z, exp - n * FLINT_BITS, rnd);
}

void
mpr_get_mpfr(mpfr_t y, const mpr_t x, mpfr_rnd_t rnd)
{
    _mpr_get_mpfr_signed(y, x->d, x->exp, x->size, x->sign, rnd);
}