#include "mpr.h"
#include "ufloat.h"

typedef struct
{
    mp_ptr d;
    long exp;
    long size;
    long alloc;
    int sign;
}
mpr_struct;

typedef struct
{
    mpr_struct mid;
    ufloat_struct rad;
}
mprb_struct;

typedef mprb_struct mprb_t[1];
typedef mprb_struct * mprb_ptr;
typedef const mprb_struct * mprb_srcptr;

#define MPRB_SIGN_PLUS 0
#define MPRB_SIGN_MINUS 1

void mprb_init(mprb_t x, long bits);
void mprb_clear(mprb_t x);

void mprb_debug(const mprb_t x);
void mprb_randtest(mprb_t x, flint_rand_t state, long emin, long emax);

long mprb_get_mid_mpz_2exp(mpz_t a, const mprb_t x);

void mprb_set_mpfr(mprb_t x, const mpfr_t v);
void mprb_get_rad_mpfr(mpfr_t r, const mprb_t x);
void mprb_get_mid_mpfr(mpfr_t x, const mprb_t v, mpfr_rnd_t rnd);
void mprb_get_interval_mpfr(mpfr_t a, mpfr_t b, const mprb_t x);
int mprb_contains_mpfr(const mprb_t x, const mpfr_t f);

/* XXX: assumes nonzero! */
static __inline__ void
_mpr_get_ufloat(ufloat_t u, mp_srcptr d, long size, long exp)
{
    /* since we truncate, we need to round up */
    u->man = (d[size - 1] >> (FLINT_BITS - UFLOAT_PREC)) + 1UL;
    u->exp = exp;

    /* adjust for carry (very unlikely!) */
    if (u->man >= (1UL << UFLOAT_PREC))
    {
        u->man = (u->man >> 1) + 1UL;
        u->exp++;
    }
}

void mprb_add(mprb_t z, const mprb_t x, const mprb_t y);
void mprb_mul(mprb_t z, const mprb_t x, const mprb_t y);
